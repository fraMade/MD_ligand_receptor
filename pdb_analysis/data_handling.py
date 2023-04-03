import pandas as pd
import numpy as np
import base64
import json
import re
import io

'''
Given the interactions data, the function count how many interaction the ligand does per bond-type.
Return a pandas dataframe with bond-types as rows and the number of timestamps associated with those. 
'''
def getBondPermanence(json_data, lig_atoms):
    bond_permanence = dict()   

    for bond_type in json_data:
        timestamps = set()
        for bond_id in json_data[bond_type]:
            if any(atom[0] in bond_id for atom in lig_atoms):
                timestamps.update( json_data[bond_type][bond_id]['ps'])
                         
        bond_permanence[bond_type] = len(timestamps) 
   
    df = pd.DataFrame.from_dict(bond_permanence, orient='index', columns=[ '#Timestamps'])
    df.index.name = 'BondType'
   
    return df.reset_index()


'''
Given the interactions data, the function count the interactions per bond-type done by each of the given atoms .
It returns a pandas dataframe with the atoms as rows and the count of interactions for each of the bond-type as columns.
'''
def get_atoms_permanence(lig_atoms:list, json_data, total_timestamp):
    atoms_permanences = {}
   
    ## iterate over the data
    for atom in lig_atoms:
        atom_id = atom[1]
        atoms_permanences[atom[1]] = {x:0 for x in json_data.keys()}
        for bond_type in json_data.keys():
            interactions = set()
            for bond_id, data in json_data[bond_type].items():
                if atom[0] in bond_id:
                    interactions.update(data['ps'])
            
            atoms_permanences[atom[1]][bond_type]= len(interactions)
    ##create the dataframe and calculate the permanences percentage                   
    df = pd.DataFrame.from_dict(atoms_permanences, orient='index', columns = json_data.keys())
    df["percentage"]= df.sum(axis=1)/(total_timestamp+1) * 100
    return df


'''
Given the interactions data, the function return a dataframe containing the list of timeframes per bond-type in which the given atom is involved at least one time .
'''
def permanence_bond(json_data, atom_id):
    
    bond_permanence = {x:set() for x in json_data.keys()}

    for bond_type in bond_permanence.keys():
        
        for bond_id in json_data[bond_type]:
             if atom_id in bond_id:
                bond_permanence[bond_type].update(json_data[bond_type][bond_id]["ps"])
       
               
    temp = []
    for k, v in bond_permanence.items():  temp+=list(zip([k]*len(v), list(map(lambda vv: int(vv),v))))          
    df = pd.DataFrame(temp, columns=['Receptor','Timestamp']).sort_values(by = 'Receptor', ascending = False)
       
    return df



'''
Given the interactions data, the function return two dataframe containing the number of interactions for each pair of (ligand-atom,aminoacid) and (ligand-atom, nuclotide).
'''
def atom_receptor2(json_data,lig_atoms, bonds_types, timestamp_num):

    nucleotid =  {x[1]:{} for x in lig_atoms}
  
    ammin = {x[1]:{} for x in lig_atoms}

    for atom in lig_atoms:  
        atom_id = atom[0]
        atom_mol = atom[1]
        considered = set()
        for bond_type in bonds_types:
            for bond,data in json_data[bond_type].items():
            

                
                if atom_id in bond:
                    
                    residue = data['restype']+data['resnr']
                    
                    for timestamp in data['ps']:

                        el = residue + str(timestamp)
                        if el not in considered:
                            if len(data['restype']) == 3:

                                if residue not in ammin[atom_mol]:
                                    ammin[atom_mol][residue] = 0

                                ammin[atom_mol][residue] += 1
                            else:
                                if residue not in nucleotid[atom_mol]:
                                    nucleotid[atom_mol][residue] = 0

                                nucleotid[atom_mol][residue] += 1 

                            considered.add(el)
        
   
    nuc = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in nucleotid.items()]),dtype=np.int32).fillna(0) 
    
    aminoacid = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in ammin.items()]), dtype=np.int32).fillna(0) 
    
    return nuc, aminoacid 



'''
Given the interactions data, the function return a dataframe containing for each given atom a list of the timeframes in which it is involved.
'''
def atom_receptor_perma(json_data,lig_atoms, timestamp_num):
    bond_permanence = {}

    for bond_type in json_data.keys():
        for bond,data in json_data[bond_type].items():

            for l in lig_atoms:

                if l[0] in bond:
                    if data['restype']+data['resnr'] not in bond_permanence:
                        bond_permanence[data['restype']+data['resnr']] = set()
                    bond_permanence[data['restype']+data['resnr']].update(filter(lambda x: int(x) <= timestamp_num*10, data['ps']))
    temp = []
    for k, v in bond_permanence.items():  temp+=list(zip([k]*len(v), list(map(lambda vv: int(vv),v))))          
    df = pd.DataFrame(temp, columns=['Receptor','Timestamp']).sort_values(by = 'Receptor', ascending = False)
    
    return df


'''
Given the interactions data, the function return a dataframe containing the list of timeframes in which the given ligand atom and receptor atom bind, per bond-type.
'''
def onelig_oneres(json_data, lig_atom, res, timestamp_num):
        
    bond_permanence = {x:set() for x in json_data.keys()}

    for bond_type in json_data.keys():
        for bond,data in json_data[bond_type].items():
            if lig_atom in bond and  data['restype']+data['resnr'] == res:
                bond_permanence[bond_type].update(filter(lambda x: int(x) <= timestamp_num*10, data['ps']))
    
    temp = []
    for k, v in bond_permanence.items():  temp+=list(zip([k]*len(v), list(map(lambda vv: int(vv),v))))          
    df = pd.DataFrame(temp, columns=['Receptor','Timestamp']).sort_values(by = 'Receptor', ascending = False)

    return df


'''
Utility function to load data.
'''
def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    
    try:
        if 'json' in filename:
            data = json.loads(decoded)

            data.pop('salt-bridges', None)
        if 'csv' in filename:
            data = pd.read_csv(io.StringIO(decoded.decode('utf-8'))).to_json(date_format='iso', orient='split')
            
    except Exception as e:
        return e

    return data
