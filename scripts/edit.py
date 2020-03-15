

def editor(pdb_file): 
    
    record = []
    with open(pdb_file, "r") as f: 
        for line in f: 
            if "HETATM" in line: 
                if "HOH" not in line: 
                    record.append(line)
                elif "ATP" in line: 
                    line = line.replace("ATP", 'LIG')
                    record.append(line)
            elif 'CONECT' not in line: 
                record.append(line)                
                
    return record
                
def printer(pdb_file, record): 
    with open(pdb_file, "w") as f: 
        for i in record: 
            f.write(i)
            
c = editor("/home/nadzhou/Downloads/1ktq.pdb")
printer("1ktq.pdb", c)
            
            
            