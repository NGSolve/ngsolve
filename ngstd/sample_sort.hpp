#ifndef SAMPLE_SORT_HPP_
#define SAMPLE_SORT_HPP_

/**************************************************************************/
/* File:   sample_sort.hpp                                                */
/* Author: Bernd Schwarzenbacher, Joachim Schoeberl                       */
/* Date:   Nov 2017                                                       */
/**************************************************************************/

#include <random>

namespace ngstd
{

  template <typename T, typename TI>
void SampleSortI(FlatArray<T> data, FlatArray<TI> index)
{
  static Timer Tsample_sort("Sample Sort");
  RegionTimer Rsample_sort(Tsample_sort);

  size_t n = index.Size();
  int nr_buckets = ngstd::TaskManager::GetNumThreads();
  int sample_size = nr_buckets-1;
  int over_sample = 10;
  int over_sample_size = over_sample * sample_size;

  Array<TI> over_sampled_ind(over_sample_size);

  std::random_device rd;
  std::mt19937 engine(rd());
  std::uniform_int_distribution<int> dist(0, index.Size()-1);
  for (auto i : Range(over_sample_size))
    over_sampled_ind[i] = index[dist(engine)];
  
  QuickSortI(data, over_sampled_ind);

  Array<TI> sample_ind(sample_size);
  for (auto i : Range(sample_size)) {
    sample_ind[i] = over_sampled_ind[i * over_sample];
  }

  Array<int> bucket_of_ind(n);
  ParallelFor(n, [&] (auto i) {
      int start = 0;
      int end = sample_size-1;
      int mid = (start+end)/2;
      while (start <= end) {
        mid = (start+end)/2;
        if (data[sample_ind[mid]] < data[i]) {
          start = mid + 1;
          continue;
        }
        else if (data[sample_ind[mid]] > data[i]) {
          end = mid - 1;
          continue;
        }
        else { break; }
      }
      if (start > end) {
        bucket_of_ind[i] = start;
      }
      else {
        bucket_of_ind[i] = mid;
      }
    });

  static Timer T3_2("Sample Sort - inverse index bucket map");
  T3_2.Start();
  ngstd::TableCreator<int> buckets_creator(nr_buckets);
  for (; !buckets_creator.Done(); buckets_creator++) {
    ParallelForRange (n, [&] (IntRange r)
                      {
                        ngstd::TableCreator<int> mycreator(nr_buckets);
                        for (; !mycreator.Done(); mycreator++)                      
                          for (auto i : r)
                            mycreator.Add (bucket_of_ind[i], i);

                        auto mytab = mycreator.MoveTable();
                        for (auto i : Range(nr_buckets))
                          buckets_creator.Add (i, mytab[i]);
                      });
  }
  
  auto table = buckets_creator.MoveTable();
  T3_2.Stop();

  static Timer T4("Sample Sort - sort buckets");
  T4.Start();
  ParallelFor(nr_buckets, [&] (auto bucket) {
    QuickSortI(data, table[bucket]);
  });
  T4.Stop();

  size_t start = 0;
  for (size_t bucket = 0; bucket < table.Size(); ++bucket) {
    size_t end = start+table[bucket].Size();
    index.Range(start, end) = table[bucket];
    start = end;
  }
}


} 

#endif  // SAMPLE_SORT_HPP_
