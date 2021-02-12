from setuptools import setup

setup(
   name='chembl_ident',
   version='1.0',
   description='ChEMBL identifier object that combines ChEMBL ID, MolRegNo and Drugbase ID. \
                Useful because no single identifier covers all compounds in ChEMBL and Drugbase. \
                Accessory functions are packaged to fill in blank identifiers when not all are known.',
   author='Tim Rozday',
   author_email='timrozday@ebi.ac.uk',
   packages=['chembl_ident'],  #same as name
)
