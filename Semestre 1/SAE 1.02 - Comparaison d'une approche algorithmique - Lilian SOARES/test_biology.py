#Fait par Lilian SOARES
import biology as bio

def test_est_base():
    assert est_base(A) == True
    assert est_base(G) == True
    assert est_base(U) == False
    print("La fonction est_base marche niquel !")

    
def test_est_adn():
    assert est_adn("ATGC") == True
    assert est_adn("AUBC") == False
    print("La fonction est_adn marche niquel !")
    
    
def test_arn():
    assert arn("ATGTCAAA") == 'AUGUCAAA'
    assert arn('ATBCG') == None
    print("La fonction arn marche niquel !")
    
    
def test_arn_to_codons():
    assert arn_to_codons("CGUUAGGGG") == ['CGU', 'UAG', 'GGG']
    assert arn_to_codons("CGUUAGGG") == ['CGU', 'UAG']
    assert arn_to_codons("") == []
    print("La fonction arn_to_codons marche niquel !")
    
    
def test_load_dico_codons_aa():
    assert bio.load_dico_codons_aa('data/codons_aa.json')=={'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine', 'UUA': 'Leucine', 'UUG': 'Leucine', 'CUU': 'Leucine','CUC': 'Leucine','CUA': 'Leucine','CUG': 'Leucine','AUU': 'Isoleucine','AUC': 'Isoleucine','AUA': 'Isoleucine','AUG': 'Methionine','GUU': 'Valine','GUC': 'Valine','GUA': 'Valine','GUG': 'Valine','UCU': 'Serine','UCC': 'Serine','UCA': 'Serine','UCG': 'Serine','CCU': 'Proline','CCC': 'Proline','CCA': 'Proline','CCG': 'Proline','ACU': 'Threonine','ACC': 'Threonine','ACA': 'Threonine','ACG': 'Threonine','GCU': 'Alanine','GCC': 'Alanine','GCA': 'Alanine','GCG': 'Alanine','UAU': 'Tyrosine','UAC': 'Tyrosine','CAU': 'Histidine','CAC': 'Histidine','CAA': 'Glutamine','CAG': 'Glutamine','AAU': 'Asparagine','AAC': 'Asparagine','AAA': 'Lysine','AAG': 'Lysine','GAU': 'Aspartic acid','GAC': 'Aspartic acid','GAA': 'Glutamic acid','GAG': 'Glutamic acid','UGU': 'Cysteine','UGC': 'Cysteine','UGG': 'Tryptophan','CGU': 'Arginine','CGC': 'Arginine','CGA': 'Arginine','CGG': 'Arginine','AGU': 'Serine','AGC': 'Serine','AGA': 'Arginine','AGG': 'Arginine','GGU': 'Glycine','GGC': 'Glycine','GGA': 'Glycine','GGG': 'Glycine'}
    assert bio.load_dico_codons_aa('data/codons_aa.json')!={}
    print("Test de la fonction load_dico_codons_aa : Ok")

    
def test_codons_stop():
    assert codons_stop(load_dico_codons_aa('data/codons_aa.json')) == ['UAA', 'UAG', 'UGA']
    print("La fonction codons_stop marche niquel !")
    
    
def test_codons_to_aa():
    assert codons_to_aa(["CGU", "AAU", "UAA", "GGG", "CGU"], load_dico_codons_aa('data/codons_aa.json')) == ['Arginine', 'Asparagine']
    assert codons_to_aa(["CGU", "AAU", "GGG", "CGU"], load_dico_codons_aa('data/codons_aa.json')) == ['Arginine', 'Asparagine', 'Glycine', 'Arginine']
    assert codons_to_aa([], load_dico_codons_aa('data/codons_aa.json')) == []
    print("La fonction codons_to_aa marche niquel !")
    
    
def test_nextIndice():
    assert nextIndice(["bonjour", "coucou", "bye", "ciao", "hello"], 0, ["hello", "bye"]) == 2
    assert nextIndice(["bonjour", "coucou", "bye", "ciao", "hello"], 0, ["bonjour"]) == 0
    assert nextIndice(["bonjour", "coucou", "bye", "ciao", "hello"], 0, []) == 5
    assert nextIndice(["bonjour", "coucou", "bye", "ciao", "hello"], 0, ["hello"]) == 4
    assert nextIndice([], 0, ["hello", "bye"]) == 0
    assert nextIndice(["bonjour", "coucou", "bye", "ciao", "hello"], 3, ["hello", "bye"]) == 4
    print("La fonction nextIndice marche niquel !")
    
    
def test_decoupe_sequence():
    assert decoupe_sequence(["val1", "début", "val2", "val3", "end", "val4", "fin", "begin", "val5", "fin", "val6", "début"], ["début", "begin"], ["fin", "end"]) == [['val2', 'val3'], ['val5']]
    assert decoupe_sequence(["val1", "début", "val2", "val3", "end", "val4", "fin", "begin", "val5", "fin", "val6", "début"], ["début", "begin"], ["fin"]) == [['val2', 'val3', 'end', 'val4'], ['val5']]
    assert decoupe_sequence(["val1", "début", "val2", "val3", "end", "val4", "fin", "begin", "val5", "fin", "val6", "début"], ["begin"], ["fin", "end"]) == [['val5']]
    assert decoupe_sequence(["val1", "début", "val2", "val3", "end", "val4", "fin", "begin", "val5", "fin", "val6", "début"], [], ["fin", "end"]) == []
    assert decoupe_sequence(["val1", "début", "val2", "val3", "end", "val4", "fin", "begin", "val5", "fin", "val6", "début"], ["début", "begin"], []) == [['val2','val3','end','val4','fin','begin','val5','fin','val6','début']]
    assert decoupe_sequence([], ["début", "begin"], ["fin", "end"]) == []
    assert decoupe_sequence([], [], []) == []
    print("La fonction decoupe_sequence marche niquel !")
    

def test_codons_to_seq_codantes():
    assert codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC", "CGU", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json"))
    assert codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json")) == [['CGU', 'AUG', 'AAU'], ['GGG', 'CCC']]
    assert codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "GGG", "CCC", "CGU", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json")) == [['CGU', 'AUG', 'AAU']]
    assert codons_to_seq_codantes(["CGU", "UUU", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC", "CGU", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json")) == [['AAU'], ['GGG', 'CCC', 'CGU']]
    assert codons_to_seq_codantes([], load_dico_codons_aa("data/codons_aa.json"))== []
    assert codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AAU", "UAA", "AUG", "GGG", "CCC", "CGU", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json")) == [['CGU', 'AAU'], ['GGG', 'CCC', 'CGU']]
    assert codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC", "CGU", "UAG", "GGG"], ["GGG"]) == []
    print("La fonction codons_to_seq_codantes marche niquel")
    
    
def test_seq_codantes_to_seq_aas():
    assert seq_codantes_to_seq_aas(codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC", "CGU", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json")), load_dico_codons_aa("data/codons_aa.json")) == [['Arginine', 'Methionine', 'Asparagine'], ['Glycine', 'Proline', 'Arginine']]
    assert seq_codantes_to_seq_aas(codons_to_seq_codantes(["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC", "GGU", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json")), load_dico_codons_aa("data/codons_aa.json")) == [['Arginine', 'Methionine', 'Asparagine'], ['Glycine', 'Proline', 'Glycine']]
    assert seq_codantes_to_seq_aas(codons_to_seq_codantes(["CGU", "UUU", "AUG", "CG", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC", "CGU", "UAG", "GGG"], load_dico_codons_aa("data/codons_aa.json")), load_dico_codons_aa("data/codons_aa.json")) == [['Methionine', 'Asparagine'], ['Glycine', 'Proline', 'Arginine']]
    assert seq_codantes_to_seq_aas(codons_to_seq_codantes([], load_dico_codons_aa("data/codons_aa.json")), load_dico_codons_aa("data/codons_aa.json")) == []
    print("La fonction codons_to_seq_aas marche niquel")
    
def test_adn_encore_molecule():
    assert adn_encode_molecule("CGTTTTATGCGTATGAATTAAATGGGGCCCCGTTAGGGG", load_dico_codons_aa("data/codons_aa.json"), ["Glycine", "Proline", "Arginine"]) == True
    assert adn_encode_molecule("", load_dico_codons_aa("data/codons_aa.json"), ["Glycine", "Proline", "Arginine"]) == False
    assert adn_encode_molecule("GG", load_dico_codons_aa("data/codons_aa.json"), ["Glycine", "Proline", "Arginine"]) == False 
    assert adn_encode_molecule("CGGG", load_dico_codons_aa("data/codons_aa.json"), ["Glycine", "Proline"]) == False
    assert adn_encode_molecule("CGTTTTATGCGTATGAATTAAATGGGGCCCCGTTAGGGG", load_dico_codons_aa("data/codons_aa.json"), []) == False
    assert adn_encode_molecule("CGTTTTATG", load_dico_codons_aa("data/codons_aa.json"), ["Arginine"]) == True
    print("La fonction adn_encore_molecule marche niquel !")