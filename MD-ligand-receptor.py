from mpi4py import MPI
import numpy as np
import subprocess
import os
import csv
import shutil
import json
import argparse
import sys
from timeit import default_timer as timer
import xml.etree.ElementTree as ET

#pdb per cycle
N_PDB = 100

def pdb_pipeline(topol, trj, i_start, i_end, restart, dnaLig, time_limit, worker_id, size, out_path):

    #set time limit in seconds for each process
    time_limit_s = time_limit * 60 * 60

    set_error = 0
     
    table = {"hbonds":dict(), "pi-stacks":dict(), "salt-bridges":dict(), "hydrophobic-interactions":dict(), "water-bridge":dict()}

    pdb_path = os.path.join(out_path, "work", "worker_{}".format(worker_id), "pdb")

    json_path = os.path.join(out_path, "work", "worker_{}".format(worker_id), "json_paths")

    save_path = os.path.join(out_path, "work", "worker_{}".format(worker_id), "save.json")

    try:
        if not os.path.exists(pdb_path):
            if restart:
                print("[ERROR]: missing pdb_path for worker_{} needed for restart".format(worker_id), file=sys.stderr)
                set_error = 1
            else:
                os.makedirs(pdb_path)

        if not os.path.exists(json_path):
            if restart:
                print("[ERROR] missing json_path for worker_{} needed for restart".format(worker_id), file=sys.stderr)
                set_error = 1
            else:
                os.makedirs(json_path)


    except OSError as os_error:
        print("[ERROR] Can't create directory {} for worker {}: {}".format(worker_id , os_error.filename, os_error.strerror), file=sys.stderr)
        set_error = 1

    
#--------------RESTART---------------#

    json_path = os.path.join(json_path, "bonds.json")

    #if restart, load the old interactions table and save points
    try:
        if restart:
            if not os.path.exists(save_path):
                print("[ERROR] missing save_path for worker_{} needed for restart".format(worker_id), file=sys.stderr)
            
                set_error = 1
                
                return table, 1

            else:

                i_start, i_end = getStatus(save_path)

                with open(json_path, 'r') as f:
                    table = json.load(f)

                print("Worker_{}: restart with start-end -> {}".format(worker_id,str((i_start, i_end))))

    except OSError as os_error:
        print("[ERROR] cant load json file {} from wroker {}: {}".format(os_error.filename, worker_id, os_error.strerror), file=sys.stderr)
        set_error = 1

#----------------------------------- #

    elapsed_time = 0
    #print("Rank: {} | times: {} - {}".format(worker_id, i_start, i_end))
    pdb_start = i_start

    #set plip comand if dna/rna is to be trated as ligand or receptor
    plip_comand = 'singularity run -H $PWD plip.simg --dnareceptor'
    if dnaLig:
        plip_comand = 'singularity run -H $PWD plip.simg'

    while(elapsed_time < time_limit_s and pdb_start <= i_end):

        if set_error == 0:

            try:

                start_timer = timer()

                pdb_end = pdb_start + N_PDB

                if pdb_end > i_end:
                    pdb_end = i_end

                subprocess.run("printf \"0\n\" | gmx -quiet trjconv -s "+topol+" -f "+trj+" -o "+pdb_path+"/temp_.pdb -sep -b "+str(pdb_start)+" -e "+str(pdb_end), universal_newlines=True, capture_output=True, shell=True, check=True)

                with os.scandir(pdb_path) as directory:
                    for entry in directory:
                        if entry.name.endswith(".pdb") and entry.is_file():

                            subprocess.run("sed -i -E \"/MODEL\ +/c\MODEL        1\" "+entry.path, universal_newlines=True, capture_output=True, check=True, shell=True)

                            with open(entry.path, 'r') as pdb_file:
                                for i, line in enumerate(pdb_file):
                                    if i == 1:
                                        ps = int(line.split()[5][:-6])
                                        break

                            plip_path = os.path.join(out_path, "work", "worker_{}".format(worker_id),"plip", entry.name[:-4] )

                            subprocess.run(plip_comand +" -f "+entry.path+" -o "+plip_path+" -x -q", universal_newlines=True, capture_output=True, shell=True, check=True)

                            os.remove(entry.path)

                            parseXML(plip_path+"/report.xml", table, ps)

                            shutil.rmtree(plip_path)

                end_time = timer()
                elapsed_time += end_time - start_timer

                pdb_start += N_PDB + 10

            except subprocess.CalledProcessError as proc_err:
                print("[ERROR] cmd error from worker {}:\n {}".format(worker_id, proc_err.stderr), file=sys.stderr)

                set_error = 1

            except OSError as os_error:
                print("[ERROR]: OSError on {} from worker {}: {}" % (os_error.filename, worker_id, os_error.strerror), file=sys.stderr)

                set_error = 1

        if set_error != 0:
            print("Error! Saving status for worker {} \nUse the restart option [-r] with the same number of processes to allow a retry and complete the analysis.".format(worker_id))
            saveStatus(save_path, pdb_start, pdb_end, worker_id, size)
            saveToJSON(table, json_path)

            return table, 1


#------------out while----------#

    #save status
    saveStatus(save_path, pdb_start, pdb_end, worker_id, size)

    #salve table to json file
    saveToJSON(table, json_path)
    
    #if beyond execution time limit then notify the user and return error code to main (to avoid the deletion of work folder)
    if elapsed_time > time_limit_s:
        if worker_id==0:
            print('Execution time limit reached ({}). Stopping execution and saving current progress!\n\nTo complete the analysis execute again the program with the same command line plus the restart option [-r].'.format(time_limit))
        
        return table, 1

    print("Rank {} | Execution Time: {}".format(worker_id, elapsed_time))

    return table, 0

def saveToCSV(dir_path, table):
    '''Save the interaction table in CSV format for all the supported interactions.'''

    try:

        for bond_type in table.keys():

            csv_path = os.path.join(dir_path, "{}.csv".format(bond_type))

            with open(csv_path, 'w') as csv_file:

                # csv_file.seek(0)
                # csv_file.truncate()

                for entry in table[bond_type]:

                    csv_writer = csv.DictWriter(csv_file, table[bond_type][entry])

                    if csv_file.tell() == 0:
                        csv_writer.writeheader()

                    csv_writer.writerow(table[bond_type][entry])

    except OSError as os_error:
        print("[ERROR] OSError while saving data in csv files {}: {}".format(os_error.filename, os_error.strerror), file=sys.stderr)

    return

def saveStatus(save_path, start, end, worker_id, size):
    '''Save status: last timeframe analyzed and the last timeframe to analyze.'''
    save_stats = {'start':int(start), 'end':int(end), 'worker_id':int(worker_id), 'size':int(size)}

    try:
        with open(save_path, 'w') as f:
            json.dump(save_stats, f)
    except OSError:
        print("[ERROR] Unable to save status on file worker {}: {}".format(worker_id, save_path), file=sys.stderr)
    return

def getStatus(save_path):
    """Load the save status: last timeframe analyzed and the last timeframe to analyze."""

    with open(save_path, 'r') as f:

        save_stats = json.load(f)

        start = save_stats['start']
        end = save_stats['end']

    return start, end

def saveToJSON(table, path):
    '''Save interactions table to JSON file'''
    try:
        with open(path, 'w') as f:
            json.dump(table, f)
    except OSError as os_error:
        print("[ERROR] JSON save has failed {} : {}".format(os_error.filename, os_error.strerror), file=sys.stderr)
    return


def hbondsFields(root, table, ps):
    '''Inserts the relevant data from hbond interactions'''

    for field in root.findall('bindingsite/interactions/hydrogen_bonds/hydrogen_bond'):

        acceptoridx = field.find('acceptoridx').text
        donoridx =   field.find('donoridx').text

        # interaction id
        bondID = str((acceptoridx, donoridx))

        # if it's a new interaction
        if bondID not in table["hbonds"].keys():
            table["hbonds"][bondID] = { "acc/donor-id": (acceptoridx, donoridx), "ps":[], "protisdon":bool(field.find('protisdon').text), "restype":field.find('restype').text, "resnr":field.find('resnr').text , "dist_h-a":[], "dist_d-a":[],
                 "don_angle":[],"donortype":field.find('donortype').text , "acceptortype":field.find('acceptortype').text}

        for x in table["hbonds"][bondID].keys():
            if isinstance(table["hbonds"][bondID][x], list) and x!="ps":
                table["hbonds"][bondID][x].append(float(field.find(x).text))

        # the interaction timestamp in picoseconds
        table["hbonds"][bondID]["ps"].append(ps)
    return

def saltbridgesFields(root, table, ps):
    '''Inserts the relevant data from salt-bridges interactions'''

    for field in root.findall('bindingsite/interactions/salt_bridges/salt_bridge'):

        prot = []
        lig = []

        for x in field.find('lig_idx_list').iter():
            if(x.text.startswith('\n')):
                continue
            lig.append(x.text)

        for x in field.find('prot_idx_list').iter():
            if(x.text.startswith('\n')):
                continue
            prot.append(x.text)

        # interaction id
        bondID = str((prot, lig))

        # if it's a new interaction
        if bondID not in table["salt-bridges"].keys():
            table["salt-bridges"][bondID] = {"lig/prot-id": (lig, prot), "ps":[], "restype":field.find('restype').text  , "resnr": field.find('resnr').text , "dist":[]}

        table["salt-bridges"][bondID]["dist"].append(float(field.find("dist").text))

        # the interaction timestamp in picoseconds
        table["salt-bridges"][bondID]["ps"].append(ps)

    return

def pistacksFields(root, table, ps):
    '''Inserts the relevant data from pi-stacks interactions'''

    for field in root.findall('bindingsite/interactions/pi_stacks/pi_stack'):

        prot = []
        lig = []

        for x in field.find('lig_idx_list').iter():
            if(x.text.startswith('\n')):
                continue
            lig.append(x.text)

        for x in field.find('prot_idx_list').iter():
            if(x.text.startswith('\n')):
                continue
            prot.append(x.text)


        # interaction id
        bondID = str((lig, prot))

        # if it's a new interaction
        if bondID not in table["pi-stacks"].keys():
            table["pi-stacks"][bondID] = {"lig/prot-id":(lig, prot), "ps":[], "restype":field.find('restype').text , "resnr": field.find('resnr').text , "centdist":[], "angle":[] }


        table["pi-stacks"][bondID]["centdist"].append(float(field.find("centdist").text))
        table["pi-stacks"][bondID]["angle"].append(float(field.find("angle").text))

        # the interaction timestamp in picoseconds
        table["pi-stacks"][bondID]["ps"].append(ps)

    return

def hydrophobicinteractionsFields(root, table, ps):
    '''Inserts the relevant data from hydrophobic interactions'''

    for field in root.findall('bindingsite/interactions/hydrophobic_interactions/hydrophobic_interaction'):

        ligcarbonidx = field.find('ligcarbonidx').text
        protcarbonidx = field.find('protcarbonidx').text

        # interaction id
        bondID = str((ligcarbonidx, protcarbonidx))

        # if it's a new interaction
        if bondID not in table["hydrophobic-interactions"].keys():
            table["hydrophobic-interactions"][bondID] = {"lig/prot-id":(ligcarbonidx, protcarbonidx), "ps":[], "restype":field.find('restype').text , "resnr": field.find('resnr').text ,
            "dist":[] }

        table["hydrophobic-interactions"][bondID]["dist"].append(float(field.find('dist').text))


        # the interaction timestamp in picoseconds
        table["hydrophobic-interactions"][bondID]["ps"].append(ps)

    return

def waterBridges(root, table, ps):
    '''Inserts the relevant data from water bridges interactions'''
   
    for field in root.findall('bindingsite/interactions/water_bridges/water_bridge'):

        donor_idx = field.find('donor_idx').text
        acceptor_idx = field.find('acceptor_idx').text
        water_idx = field.find('water_idx').text

        # interaction id
        bondID = str((acceptor_idx, donor_idx, water_idx))

        # if it's a new interaction
        if bondID not in table["water-bridge"].keys():
            table["water-bridge"][bondID] = {"acc/don/wat-id":(acceptor_idx, donor_idx, water_idx), "ps":[], "restype":field.find('restype').text , "resnr":field.find('resnr').text,
            "dist_a-w":[], "dist_d-w":[], "don_angle": [], "water_angle":[], "donortype":field.find('donortype').text, "acceptortype":field.find('acceptortype').text }

        for x in table["water-bridge"][bondID].keys() :
                if isinstance(table["water-bridge"][bondID][x], list) and x != "ps":
                    table["water-bridge"][bondID][x].append(float(field.find(x).text))

        # the interaction timestamp in picoseconds
        table["water-bridge"][bondID]["ps"].append(ps)

    return



def parseXML(xmlfile, table, ps):
    '''Fills the table with all the 4 interaction types'''
    
    # create element tree object
    tree = ET.parse(xmlfile)

    # get root element
    root = tree.getroot()


    hbondsFields(root, table, ps)

    hydrophobicinteractionsFields(root, table, ps)

    saltbridgesFields(root, table, ps)

    pistacksFields(root, table, ps)

    waterBridges(root, table, ps)

    return

'''
table["hbonds"][bondID] = { "acc/donor-id": bondID,
                            "ps":[],
                            "protisdon": True,
                            "restype": string,
                            "resnr": string,
                            "dist_h-a":[], "dist_d-a":[],
                            "don_angle":[],
                            "donortype":[] ,
                            "acceptortype":[]}


table["salt-bridges"][bondID] = {"prot/lig-id": bondID,
                                "ps":[],
                                "restype":field.find('restype').text  ,
                                "resnr": field.find('resnr').text ,
                                "dist":[]}


table["pi-stacks"][bondID] = {"prot/lig-id":bondID,
                              "ps":[],
                              "restype":field.find('restype').text ,
                              "resnr": field.find('resnr').text ,
                              "centdist":[],
                              "angle":[] }

table["hydrophobic-interactions"][bondID] = {"lig/prot-id":bondID,
                                            "ps":[],
                                            "restype":field.find('restype').text ,
                                            "resnr": field.find('resnr').text ,
                                            "dist":[] }


table["water-bridge"][bondID] = {"acc/don/wat-id":bondID,
                                "ps":[],
                                "restype":field.find('restype').text ,
                                "resnr":field.find('resnr').text,
                                "dist_a-w":[],
                                "dist_d-w":[],
                                "don_angle": [],
                                "water_angle":[],
                                "donortype":[],
                                "acceptortype":[] }
'''
def extract_atoms(topol, trj, timestamp, out_path, ligand_id, dnaLig):

    pdb_path = os.path.join(out_path, "work", "worker_0", "pdb")
    csv_path = os.path.join(out_path, "result", "csv_files")
    
    try:
        subprocess.run("printf \"0\n\" | gmx -quiet trjconv -s "+topol+" -f "+trj+" -o "+pdb_path+"/temp_.pdb -sep -b 0 -e 0", universal_newlines=True, capture_output=True, shell=True, check=True)
        file_path = os.path.join(pdb_path, 'temp_0.pdb')
        subprocess.run("sed -i -E \"/MODEL\ +/c\MODEL        1\" "+file_path, universal_newlines=True, capture_output=True, check=True, shell=True)
        
        #set the use of DNA/RNA as either ligand or receptor
        plip_comand = 'singularity run -H $PWD plip.simg --nofixfile --dnareceptor'
        if dnaLig:
            plip_comand = 'singularity run -H $PWD plip.simg --nofixfile'

        subprocess.run(plip_comand +" -f "+file_path+" -o "+pdb_path+" -q", universal_newlines=True, capture_output=True, shell=True, check=True)
        
        os.remove(file_path)
        plip_path = os.path.join(pdb_path, 'temp_0_protonated.pdb')

        protein = []
        ligand = []
        with open(plip_path, 'r') as pdb_file:
             for line in pdb_file.readlines():
                splitted_line = line.split()
                if splitted_line[0] == 'ATOM' and splitted_line[2]!= 'H':
                    if splitted_line[3] == ligand_id:
                        ligand.append(splitted_line[1:])
                    elif splitted_line[3] != 'SOL' and splitted_line[3] != 'NA':
                        protein.append(splitted_line[1:])
        
        os.remove(plip_path)
        
        with open(csv_path + "/LIG_ATOMS.csv", 'w') as csv_file:
            file_writer = csv.writer(csv_file)
            file_writer.writerow(["ATOM_ID","ATOM","RES_TYPE","CHAIN_ID","RES_NUM","X","Y","Z","OCCUP","TEMP","ELEMENT"])
            for atom in ligand:
                file_writer.writerow(atom)

        with open(csv_path + "/RCPT_ATOMS.csv", 'w') as csv_file:
            file_writer = csv.writer(csv_file)
            file_writer.writerow(["ATOM_ID","ATOM","RES_TYPE","CHAIN_ID","RES_NUM","X","Y","Z","OCCUP","TEMP","ELEMENT"])
            for atom in protein:
                file_writer.writerow(atom)

    except subprocess.CalledProcessError as proc_err:
        print("[ERROR] cmd error while extracting atoms_id: {}".format(proc_err.stderr), file=sys.stderr)
    
    except OSError as os_error:
        print("[ERROR] OSError while extracting atoms data {}: {}".format(os_error.filename, os_error.strerror), file=sys.stderr)

    return

def mergeTables(interactions_tables, json_path):
    '''Merge the workers's dictionaries in the main process.'''

    table = {"hbonds":dict(), "pi-stacks":dict(), "salt-bridges":dict(), "hydrophobic-interactions":dict(), "water-bridge":dict()}

    for dictionary in interactions_tables:

        for bond_type, bonds in dictionary.items():

            #dict inside bond_type

            for james_bond, bond_dict in bonds.items():

                #dict of single bond

                if james_bond not in table[bond_type].keys():
                    table[bond_type][james_bond] = dict()

                for field_names, field_values in bond_dict.items():
                    
                    #print(str(worker_id)+'\n')

                    #merge fields values
                    if type(field_values) == list:
                        if field_names not in table[bond_type][james_bond].keys():
                            table[bond_type][james_bond][field_names] = []

                        table[bond_type][james_bond][field_names].extend(field_values)

                    else:
                        if field_names not in table[bond_type][james_bond].keys():
                            table[bond_type][james_bond][field_names] = field_values

    saveToJSON(table, json_path)

    return table


'''
Main process assigns workload to workers processes then it checks the wait for all the workers to join.
Time limit (24h) is checked before merging the json/csv files to avoid unexpected termination.
'''
def main(topol, trj, start, end, restart, out_path, ligand_id, dnaLig, timeLimit):

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    err_code = 0
    
    work = np.arange(start, end + 10, 10, dtype=np.uintc)

    recv_size = int(((end - start) / 10) // size + 1)

    if size > work.shape[0]:
        if rank == 0:
            if(end - start <= 0):
                print("[ERROR] begin time is less or equal to end time!", file=sys.stderr)
            print("[ERROR] the program was started with an exceeding number of processes: {} for {} timestamps. Aborting execution!".format(size, work.shape[0]), file=sys.stderr)
        return 1

    recvbuff = np.ones(recv_size, dtype=np.uintc)
        
    comm.Scatterv(work, recvbuff, root=0)
   
    pdb_start =  recvbuff[0]
    pdb_end = recvbuff[-1] if recvbuff[-1] != 1 else np.max(recvbuff)

    if rank == 0:

        if(end - start <= 0):
            print("[ERROR] begin time is less or equal to end time!", file=sys.stderr)
            err_code = 1

        if not (end % 10 and start % 10) == 0:
           print("[ERROR] invalid value for begin and end time, it must be divisible by 10!", file=sys.stderr)
           err_code = 1
            
        for in_file in [topol, trj]:
            if not os.path.exists(in_file):
                print("[ERROR] input file: {} does not exist".format(in_file), file=sys.stderr)
                err_code = 1

    if not os.path.exists(out_path):
        if rank == 0:
            print('Output path: {} not found, will revert to the working directory'.format(out_path))
        out_path = '.'

    err_code = comm.bcast(err_code, root=0)

    #some error in the input check portion of the code above
    if err_code == 1:
        return 1

    table, err_code = pdb_pipeline(topol, trj, pdb_start, pdb_end, restart, dnaLig, timeLimit, rank, size, out_path)

    # in case of error interrupt execution and close the program
    if err_code == 0:
        interactions_tables = comm.gather(table, root=0)
    else:
        interactions_tables = comm.gather(-1, root=0)
        return 

    if rank == 0:
        
        csv_path = os.path.join(out_path, "result/csv_files")
        
        if not os.path.exists(csv_path):
            os.makedirs(csv_path)

        table = mergeTables(interactions_tables, os.path.join(out_path,"result/bonds_complete.json"))

        saveToCSV(csv_path, table)

        extract_atoms(topol, trj, pdb_start, out_path, ligand_id, dnaLig)

        #remove work directory if program closed succesfully
        shutil.rmtree(os.path.join(out_path,"work"))

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='MD-Ligand-Receptor, an efficient high-performance structural bioinformatics pipeline for the analysis of ligand-receptor binding interactions as a function of time, obtained from the analysis of molecular dynamics trajectory.')

    parser.add_argument('-s', help='file with the starting structure of the simulation, the molecular topology and all the simulation parameters: tpr', type=str, required=True)
    parser.add_argument('-f', help='Trajectory of the simulation: xtc', type=str, required=True)
    parser.add_argument('-b', help='Set begin time (ps)', type=int, required=True)
    parser.add_argument('-e', help='Set end time (ps)', type=int, required=True)
    parser.add_argument('-l', help='Ligand identifier in pdb file (e.g. \"LIG\")', type=str, required=True)
    parser.add_argument('-o', help='Output path for result and work directiories', default='.', type=str, required=False)
    parser.add_argument('-r', help='Restart from previous execution', action='store_true', required=False)
    parser.add_argument('--dnaligand', help='Set DNA/RNA molecules as ligands', action='store_true', required=False)
    parser.add_argument('-t', '--time', help='Set execution time limit (unit time hours)', default=23, type=int, required=False)

    args = parser.parse_args()

    main(args.s, args.f, args.b, args.e, args.r, args.o, args.l, args.dnaligand, args.time)