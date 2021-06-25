import json
import pickle
from collections import defaultdict
import cx_Oracle
from tqdm.auto import tqdm

class ChemblIdent():
    def __init__(self, chembl_id=None, drugbase_id=None, molregno=None):
        self.drugbase_id = drugbase_id
        self.molregno = molregno
        self.chembl_id = chembl_id
        
    def __repr__(self):
        return repr({'chembl_id': self.chembl_id, 'drugbase_id': self.drugbase_id, 'molregno': self.molregno})
    
    def __iter__(self):
        for k,v in self.__dict__.items():
            yield k,v
            
    def __hash__(self):
        return hash((self.chembl_id, self.drugbase_id, self.molregno))
    
    def __eq__(self, other):
        return (self.chembl_id, self.drugbase_id, self.molregno) == (other.chembl_id, other.drugbase_id, other.molregno)
    
    def __ne__(self, other):
        return not(self == other)
    
    def __lt__(self, other):
        def get_tuple(obj):
            if obj.chembl_id is None:
                chembl_id = -1
            else:
                chembl_id = int(obj.chembl_id[6:])  # chop off the "CHEMBL" bit of the ChEMBL ID
                
            if obj.drugbase_id is None:
                drugbase_id = -1
            else:
                drugbase_id = obj.drugbase_id
                
            if obj.molregno is None:
                molregno = -1
            else:
                molregno = obj.molregno
                
            return (chembl_id, drugbase_id, molregno)
        
        return get_tuple(self) < get_tuple(other)
    
    def __tuple__(self):
        return (self.chembl_id, self.drugbase_id, self.molregno)

class ChemblIndexes():
    def __init__(self, data_dir='.'):
        
        self.data_dir = data_dir
        
        try:
            self.load_indexes()
                
        except Exception as e:
            self.drugbase_id2mrn = None
            self.mrn2drugbase_id = None
            self.chembl_id2mrn = None
            self.mrn2chembl_id = None
            self.compound_parents = None
            self.compound_children = None
        
    def connect_to_chembl(self, username, password, url):
        self.chembl_db = cx_Oracle.connect(username, password, url)
#         print(chembl_db.dsn, "running version", chembl_db.version)
        
        return self.chembl_db
        
    def gen_indexes(self):
        chembl_cursor = self.chembl_db.cursor()
        
        sql_query = (   
                        "select MD.ID, MD.MOLREGNO "
                        "from DRUGBASE.MOLECULE_DICTIONARY MD "
                        "where MD.DELETED = 0 "
                    )
        chembl_cursor.execute(sql_query)
        
        self.drugbase_id2mrn = {}
        self.mrn2drugbase_id = {}
        for drugbase_id, molregno in tqdm(chembl_cursor, desc='Drugbase molecule dictionary'):
            if (drugbase_id is None) or (molregno is None):
                continue
            self.drugbase_id2mrn[drugbase_id] = molregno
            self.mrn2drugbase_id[molregno] = drugbase_id
            
        sql_query = (   
                        "select MD.CHEMBL_ID, MD.MOLREGNO "
                        "from CHEMBL.MOLECULE_DICTIONARY MD "
                    )
        chembl_cursor.execute(sql_query)
        
        self.chembl_id2mrn = {}
        self.mrn2chembl_id = {}
        for chembl_id, molregno in tqdm(chembl_cursor, desc='ChEMBL molecule dictionary'):
            if (chembl_id is None) or (molregno is None):
                continue
            self.chembl_id2mrn[chembl_id] = molregno
            self.mrn2chembl_id[molregno] = chembl_id
        
        # fetch phase data
        
        sql_query = (   
                        "select MD.MOLREGNO, MD.MAX_PHASE "
                        "from CHEMBL.MOLECULE_DICTIONARY MD "
                    )
        chembl_cursor.execute(sql_query)

        self.mrn2phase = {}
        for molregno, phase in tqdm(chembl_cursor, desc='ChEMBL trial phase'):
            if (molregno is None):
                continue
            self.mrn2phase[molregno] = phase
            
        sql_query = (   
                        "select MD.ID, MD.HIGHEST_PHASE "
                        "from DRUGBASE.MOLECULE_DICTIONARY MD "
                    )
        chembl_cursor.execute(sql_query)

        self.drugbase_id2phase = {}
        for drugbase_id, phase in tqdm(chembl_cursor, desc='Drugbase trial phase'):
            if (drugbase_id is None):
                continue
            self.drugbase_id2phase[drugbase_id] = phase
        
        # fetch drugbase sources
        
        sql_query = (   
                        "select MS.ID, MS.SOURCE "
                        "from DRUGBASE.MOLECULE_SOURCE MS "
                    )
        chembl_cursor.execute(sql_query)

        self.source_map = {}
        for source_id, source in tqdm(chembl_cursor, desc='Drugbase source map'):
            self.source_map[source_id] = source
            
        sql_query = (   
                        "select MSM.MOLECULE_DICTIONARY_ID, MSM.MOLECULE_SOURCE_ID "
                        "from DRUGBASE.MOLECULE_SOURCE_MAPPING MSM "
                        "where MSM.REMOVED = 0 "
                    )
        chembl_cursor.execute(sql_query)

        self.drugbase_id2source_id = defaultdict(set)
        for drugbase_id, source_id in tqdm(chembl_cursor, desc='Drugbase source IDs'):
            self.drugbase_id2source_id[drugbase_id].add(source_id)

        self.drugbase_id2source_id = dict(self.drugbase_id2source_id)
        
        # fetch compound hierarchy
        
        self.compound_parents = defaultdict(set)
        for child, parent in tqdm(chembl_cursor.execute("select MOLREGNO, PARENT_MOLREGNO from CHEMBL.MOLECULE_HIERARCHY"), leave=True, position=0, desc='ChEMBL hierarchy'):
            child = self.get_chembl_ident(molregno=child)
            parent = self.get_chembl_ident(molregno=parent)
            self.compound_parents[child].add(parent)
            
        for child, parent, active in tqdm(chembl_cursor.execute("select MOLECULE_DICTIONARY_ID, PARENT_MOLECULE_DICTIONARY_ID, ACTIVE_MOLECULE_DICTIONARY_ID from DRUGBASE.MOLECULE_HIERARCHY"), leave=True, position=0, desc='Drugbase hierarchy'):
            child = self.get_chembl_ident(drugbase_id=child)
            parent = self.get_chembl_ident(drugbase_id=parent)
            active = self.get_chembl_ident(drugbase_id=active)
            self.compound_parents[child].update({parent, active})
        
        self.compound_children = defaultdict(set)
        for k,vs in self.compound_parents.items():
            for v in vs:
                self.compound_children[v].add(k)
        
        self.compound_parents = dict(self.compound_parents)
        self.compound_children = dict(self.compound_children)
        
        chembl_cursor.close()
    
    def save_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
            
        with open(f"{data_dir}/drugbase_id2mrn.json", 'wt') as f:
            json.dump(self.drugbase_id2mrn, f)
        with open(f"{data_dir}/mrn2drugbase_id.json", 'wt') as f:
            json.dump(self.mrn2drugbase_id, f)
        with open(f"{data_dir}/chembl_id2mrn.json", 'wt') as f:
            json.dump(self.chembl_id2mrn, f)
        with open(f"{data_dir}/mrn2chembl_id.json", 'wt') as f:
            json.dump(self.mrn2chembl_id, f)
        with open(f"{data_dir}/compound_parents.pkl", 'wb') as f:
            pickle.dump({k.__tuple__():{v.__tuple__() for v in vs} for k,vs in self.compound_parents.items()}, f)
        with open(f"{data_dir}/compound_children.pkl", 'wb') as f:
            pickle.dump({k.__tuple__():{v.__tuple__() for v in vs} for k,vs in self.compound_children.items()}, f)
        with open(f"{data_dir}/source_map.json", 'wt') as f:
            json.dump(self.source_map, f)
        with open(f"{data_dir}/drugbase_id2source_id.json", 'wt') as f:
            json.dump({k:list(vs) for k,vs in self.drugbase_id2source_id.items()}, f)
        with open(f"{data_dir}/mrn2phase.json", 'wt') as f:
            json.dump(self.mrn2phase, f)
        with open(f"{data_dir}/drugbase_id2phase.json", 'wt') as f:
            json.dump(self.drugbase_id2phase, f)
                
    def load_indexes(self, data_dir=None):
        if data_dir is None:
            data_dir = self.data_dir
        
        with open(f"{data_dir}/drugbase_id2mrn.json", 'rt') as f:
            self.drugbase_id2mrn = {int(k):v for k,v in json.load(f).items()}
        with open(f"{data_dir}/mrn2drugbase_id.json", 'rt') as f:
            self.mrn2drugbase_id = {int(k):v for k,v in json.load(f).items()}
        with open(f"{data_dir}/chembl_id2mrn.json", 'rt') as f:
            self.chembl_id2mrn = {k:int(v) for k,v in json.load(f).items()}
        with open(f"{data_dir}/mrn2chembl_id.json", 'rt') as f:
            self.mrn2chembl_id = {int(k):v for k,v in json.load(f).items()}
        with open(f"{data_dir}/compound_parents.pkl", 'rb') as f:
            self.compound_parents = {ChemblIdent(*k):{ChemblIdent(*v) for v in vs} for k,vs in pickle.load(f).items()}
        with open(f"{data_dir}/compound_children.pkl", 'rb') as f:
            self.compound_children = {ChemblIdent(*k):{ChemblIdent(*v) for v in vs} for k,vs in pickle.load(f).items()}
        with open(f"{data_dir}/source_map.json", 'rt') as f:
            self.source_map = {int(k): v for k,v in json.load(f).items()}
        with open(f"{data_dir}/drugbase_id2source_id.json", 'rt') as f:
            self.drugbase_id2source_id = {int(k):set(vs) for k,vs in json.load(f).items()}
        with open(f"{data_dir}/mrn2phase.json", 'rt') as f:
            self.mrn2phase = {int(k):v for k,v in json.load(f).items()}
        with open(f"{data_dir}/drugbase_id2phase.json", 'rt') as f:
            self.drugbase_id2phase = {int(k):v for k,v in json.load(f).items()}
    
    def get_parents(self, obj=None, drugbase_id=None, molregno=None, chembl_id=None):
        if obj is None:
            obj = self.get_chembl_ident(drugbase_id=drugbase_id, molregno=molregno, chembl_id=chembl_id)
        if obj:
            if obj in self.compound_parents:
                return self.compound_parents[obj]
        
    def get_children(self, obj=None, drugbase_id=None, molregno=None, chembl_id=None):
        if obj is None:
            obj = self.get_chembl_ident(drugbase_id=drugbase_id, molregno=molregno, chembl_id=chembl_id)
        if obj:
            if obj in self.compound_children:
                return self.compound_children[obj]
            
    def get_molregno(self, drugbase_id=None, chembl_id=None):
        if drugbase_id in self.drugbase_id2mrn:
            return self.drugbase_id2mrn[drugbase_id]
        if chembl_id in self.chembl_id2mrn:
            return self.chembl_id2mrn[chembl_id]
        
    def get_drugbase_id(self, molregno=None):
        if molregno in self.mrn2drugbase_id:
            return self.mrn2drugbase_id[molregno]
        
    def get_chembl_id(self, molregno=None):
        if molregno in self.mrn2chembl_id:
            return self.mrn2chembl_id[molregno]
    
    def get_sources(self, drugbase_id):
        if drugbase_id in self.drugbase_id2source_id: 
            return {self.source_map[s] for s in self.drugbase_id2source_id[drugbase_id]}
        else:
            return set()
        
    def get_phase(self, obj=None, drugbase_id=None, molregno=None, chembl_id=None):
        if obj is None:
            obj = self.get_chembl_ident(drugbase_id=drugbase_id, molregno=molregno, chembl_id=chembl_id)
        if obj:
            phase1 = self.mrn2phase[obj.molregno] if obj.molregno in self.mrn2phase else None
            phase2 = self.drugbase_id2phase[obj.drugbase_id] if obj.drugbase_id in self.drugbase_id2phase else None
            max_phase = max([
                -5 if phase1 is None else phase1, 
                -5 if phase2 is None else phase2, 
            ])
            if max_phase > -5:
                return max_phase
    
    def get_chembl_ident(self, molregno=None, drugbase_id=None, chembl_id=None):
        if molregno is None:
            molregno = self.get_molregno(drugbase_id=drugbase_id, chembl_id=chembl_id)
        if chembl_id is None:
            chembl_id = self.get_chembl_id(molregno=molregno)
        if drugbase_id is None:
            drugbase_id = self.get_drugbase_id(molregno=molregno)
        
        if any([molregno,drugbase_id,chembl_id]):
            return ChemblIdent(molregno=molregno, drugbase_id=drugbase_id, chembl_id=chembl_id)
    

            
