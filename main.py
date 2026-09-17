
from pathlib import Path

import anndata as ad
import h5py

input_file = Path(
    r"E:\\PATCS\\input\\GSE246613"
    r"\\GSE246613_PembroRT_immune_R100_final.h5ad"
)

output_file = input_file.with_name(
    "GSE246613_PembroRT_immune_R100_final_updated.h5ad"
)

if not input_file.exists():
    raise FileNotFoundError(f" FileNotFound : {input_file}")


try:
    with h5py.File(input_file, mode="r") as h5:
        print("\nHDF5 : ")
        print(dict(h5.attrs))

        print("\nHDF5 : ")
        print(list(h5.keys()))

except OSError as error:
    raise RuntimeError(
        "This file is not a valid HDF5 file"
    ) from error

print("\nReading old h5ad ...")

adata = ad.read_h5ad(
    filename=input_file
)

print("\nSuccess")
print("Cells: ", adata.n_obs)
print("Genes: ", adata.n_vars)
print("X type: ", type(adata.X))
print("layers: ", list(adata.layers.keys()))
print("obsm: ", list(adata.obsm.keys()))
print("varm: ", list(adata.varm.keys()))
print("obsp: ", list(adata.obsp.keys()))
print("uns: ", list(adata.uns.keys()))
print("Is there raw: ", adata.raw is not None)

adata.obs_names = adata.obs_names.astype(str)
adata.var_names = adata.var_names.astype(str)

if not adata.obs_names.is_unique:
    print("obs_names_make_unique()")
    adata.obs_names_make_unique()

if not adata.var_names.is_unique:
    print("var_names_make_unique()")
    adata.var_names_make_unique()

print("\nNew version of h5ad ...")

adata.write_h5ad(filename=output_file,compression="gzip")

print("\nSuccess: ")
print(output_file)
print("File size : ", round(output_file.stat().st_size / 1024**3, 3), "GB")

print("\nChecking ...")

adata_check = ad.read_h5ad(
    filename=output_file,
    backed="r"
)

print("Success ")
print("Cells : ", adata_check.n_obs)
print("Genes : ", adata_check.n_vars)

adata_check.file.close()

