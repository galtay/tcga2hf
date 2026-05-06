"""Build pipeline that downloads public TCGA data from the NCI GDC and
publishes it as the `gabrielaltay/tcga-{patients,tabular}-open` HuggingFace
datasets. The user-facing layer is the `tcga2hf-pipeline` CLI; for the
read-side typed models / pyarrow schemas, install the companion `tcga2hf`
package.
"""

__version__ = "0.1.0"
