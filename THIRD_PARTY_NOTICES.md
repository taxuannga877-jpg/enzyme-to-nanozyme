# Third-party web assets

The optional Flask workbench bundles browser assets so it can run locally
without a CDN. Versions, distributed license texts, and SHA-256 hashes were
frozen from the upstream package releases listed below.

| Asset | Version | License | Bundled file SHA-256 |
| --- | --- | --- | --- |
| [Bootstrap](https://www.npmjs.com/package/bootstrap/v/4.3.1) | 4.3.1 | MIT | CSS `f9dcd8538e4a5258dc86f0fbd965802ad09ee4afef515f4fe1053e4c4b19f6bf`; JS `0a34a87842c539c1f4feec56bba982fd596b73500046a6e6fe38a22260c6577b` |
| [jQuery slim](https://www.npmjs.com/package/jquery/v/3.3.1) | 3.3.1 | MIT | `dde76b9b2b90d30eb97fc81f06caa8c338c97b688cea7d2729c88f529f32fbb1` |
| [Popper.js](https://www.npmjs.com/package/popper.js/v/1.14.7) | 1.14.7 | MIT | `66f3a07e1fa9b64a686b66381e4458dbc8abf3dbbff954720c4eec07b84411c2` |
| [Font Awesome Free](https://www.npmjs.com/package/@fortawesome/fontawesome-free/v/6.0.0-beta3) | 6.0.0-beta3 | Icons CC BY 4.0; fonts SIL OFL 1.1; code MIT | CSS `a361e7885c36bacb3fd9cb068da207c3b9329962cac022d06e28923939f575e8`; individual font hashes are recorded by `BUILD_PROVENANCE.json` |
| [3Dmol.js](https://www.npmjs.com/package/3dmol/v/2.5.5) | 2.5.5 | BSD-3-Clause with incorporated GLmol, Three.js, and jQuery notices | JS `f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016`; webpack notice `ae3bfc688d0c9687b76e0ecc7fece0b393bcb4acf2aee09e19dda697eaa10b16` |

## Distributed license files

The complete upstream texts carried with this repository are:

- `enzyme_viewer/static/js/3Dmol-min.js.LICENSE.txt`;
- `enzyme_viewer/static/vendor/licenses/3Dmol-2.5.5-BSD-3-Clause.txt`;
- `enzyme_viewer/static/vendor/licenses/Bootstrap-4.3.1-MIT.txt`;
- `enzyme_viewer/static/vendor/licenses/jQuery-3.3.1-MIT.txt`;
- `enzyme_viewer/static/vendor/licenses/Popper.js-1.14.7-MIT.txt`; and
- `enzyme_viewer/static/vendor/licenses/Font-Awesome-Free-6.0.0-beta3.txt`.

The minified-file copyright headers must not be removed. Brand icons remain
trademarks of their respective owners; inclusion does not imply endorsement.
These third-party notices do not determine the still-unresolved license for
E2N's original code, data, figures, or documentation.
