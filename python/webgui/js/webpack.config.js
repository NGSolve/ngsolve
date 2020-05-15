const path = require('path');
const version = require('./package.json').version;

// Custom webpack rules
const rules = [
  { test: /\.ts$/, loader: 'ts-loader' },
  { test: /\.js$/, loader: 'source-map-loader' },
  { test: /\.css$/, use: ['style-loader', 'css-loader']}
];

// Packages that shouldn't be bundled but loaded at runtime
const externals = ['@jupyter-widgets/base'];
const externals1 = {
'three': "https://cdn.jsdelivr.net/npm/three@0.115.0/build/three.min.js",
'dat.gui': "https://cdnjs.cloudflare.com/ajax/libs/dat-gui/0.7.7/dat.gui.js"
};


const resolve = {
  // Add '.ts' and '.tsx' as resolvable extensions.
  extensions: [".webpack.js", ".web.js", ".ts", ".js"]
};

module.exports = function(env, argv)
{
  const build_dir = path.resolve(env.outdir);
  return [
  /**
   * Notebook extension
   *
   * This bundle only contains the part of the JavaScript that is run on load of
   * the notebook.
   */
  {
    target: 'web',
    entry: './src/extension.ts',
    output: {
      filename: 'index.js',
      path: path.resolve(build_dir, 'nbextension', 'static'),
      libraryTarget: 'amd'
    },
    module: {
      rules: rules
    },
    devtool: 'source-map',
    externals,
    resolve,
  },

  /**
   * Embeddable ngoslve_jupyter_widgets bundle
   *
   * This bundle is almost identical to the notebook extension bundle. The only
   * difference is in the configuration of the webpack public path for the
   * static assets.
   *
   * The target bundle is always `dist/index.js`, which is the path required by
   * the custom widget embedder.
   */
  {
    target: 'web',
    entry: './src/index.ts',
    output: {
      filename: 'index.js',
      path: path.resolve(build_dir, 'dist'),
      libraryTarget: 'amd',
      library: "ngsolve_jupyter_widgets",
      publicPath: 'https://unpkg.com/ngsolve_jupyter_widgets@' + version + '/dist/'
    },
    devtool: 'source-map',
    module: {
        rules: rules
    },
    externals,
    resolve,
  },


  /**
   * Documentation widget bundle
   *
   * This bundle is used to embed widgets in the package documentation.
   */
  {
    target: 'web',
    entry: './src/index.ts',
    output: {
      filename: 'embed-bundle.js',
      path: path.resolve(build_dir, 'docs', 'source', '_static'),
      library: "ngsolve_jupyter_widgets",
      libraryTarget: 'amd'
    },
    module: {
      rules: rules
    },
    devtool: 'source-map',
    externals,
    resolve,
  },

  /**
   * Standalone html bundle
   *
   * This bundle is used to create html output files directly from python console
   */
  {
    target: 'web',
    entry: './src/scene.ts',
    output: {
      filename: 'standalone.js',
      path: path.resolve(build_dir),
      library: "ngsolve_jupyter_widgets",
      libraryTarget: 'amd'
    },
    module: {
      rules: rules
    },
    devtool: 'source-map',
    externals: externals,
    resolve,
  }

] };
