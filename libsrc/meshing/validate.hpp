#ifndef VALIDATE_HPP
#define VALIDATE_HPP



void GetPureBadness(Mesh & mesh, Array<double> & pure_badness,
		    const BitArray & isnewpoint);
double Validate(const Mesh & mesh, Array<ElementIndex> & bad_elements,
		const Array<double> & pure_badness, 
		double max_worsening, const bool uselocalworsening,
		Array<double> * quality_loss = NULL);
void RepairBisection(Mesh & mesh, Array<ElementIndex> & bad_elements, const BitArray & isnewpoint, Refinement & refinement,
		     const Array<double> & pure_badness, 
		     double max_worsening, const bool uselocalworsening,
		     const Array< Array<int,PointIndex::BASE>* > & idmaps);

#endif // VALIDATE_HPP
