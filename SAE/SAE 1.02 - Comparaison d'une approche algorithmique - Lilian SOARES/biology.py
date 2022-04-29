def est_base(c):
    if c == 'A' or c == 'T' or c == 'G' or c == 'C':
        return True
    else:
        return False


def est_adn(s):
    i = 0
    while i < len(s) :
        c = s[i]
        if est_base(c) == True:
            i = i + 1
        elif est_base(c) == False:
            return False
            
    return True


def arn(adn):
    s = adn
    if est_adn(s) == False :
        print("Le séquencage ADN n'est pas bon")
        return None
        
    i = 0
    arn = ""
    while i < len(adn):
        if adn[i] == "T":
            arn = arn + "U"
        else :
            arn = arn + adn [i]
        i = i + 1
        
    return arn


def arn_to_codons(arn):
    i = 0
    codons = []
    mp_3 = len(arn) - len(arn)%3
    
    while i < mp_3 :
        codons.append(arn[i] + arn[i+1] + arn[i+2])
        i = i + 3
    return codons


def load_dico_codons_aa(filename):
    fichier = open(filename,"r")
    codons_json = fichier.read()
    dico_adn = loads(codons_json)
    fichier.close()
    return dico_adn


def codons_stop(dico):
    stop = []
    bases = "AUGC"
    i = 0
    while i < 4:
        j = 0
        while j < 4:
            k = 0
            while k < 4:
                if bases[i] + bases[j] + bases[k] not in dico:
                    stop.append(bases[i] + bases[j] + bases[k])
                k += 1
            j += 1
        i += 1
        
    return stop


def codons_to_aa(codons, dico):
    i = 0
    j = 0
    triplés = []
    acides = []
    while i < len(codons):
        if codons[i] not in dico :
            i = len(codons)
        else :
            triplés.append(codons[i])
        i = i + 1

    clés = list(dico)
    i = 0
    while i < len(triplés):
        while j < len(clés):
            if triplés[i] == clés[j]:
                acides.append(list(dico.values())[j])
                j = len(clés)
            else :
                j = j + 1

        i = i + 1
        j = 0

    return acides


def nextIndice(tab, ind, elements) :
    while ind<len(tab) :
        i = 0
        while i < len(elements) :
            if tab[ind] == elements[i] :
                return ind
            i = i + 1
        ind = ind + 1
    return ind


def decoupe_sequence(seq, start, stop):
    i = 0
    k = 0
    z = 0
    compteur = 0
    for g in range(len(seq)):
        if seq[g] in start:
            compteur = compteur + 1
          
    liste = [[]for j in range(compteur)]
    while i < len(start)+1 :
        a = nextIndice(seq, k, start)
        b = nextIndice(seq, k, stop)
        if a > b :
            k = a
            b = nextIndice(seq, k, stop)
        k = a + 1
        
        while k < b :
            liste[z].append(seq[k])              
            k = k + 1
        z = z + 1    
        k = b    
        i = i + 1    
        
    i = 0
    while i < len(liste):
        if liste[i] == []:
            liste.pop(i)
        else : 
            i = i + 1
            
    return liste


def codons_to_seq_codantes(seq, dico) :
    stop = codons_stop(dico)
    seq_codantes = decoupe_sequence(seq, ["AUG"], stop)
    return seq_codantes


def seq_codantes_to_seq_aas(seq_codantes, dico):
    i = 0
    j = 0
    k = 0
    clés = list(dico)
    tab = [[]for g in range(len(seq_codantes))]
    
    while i < len(seq_codantes):
        while j < len(seq_codantes[i]):
            while k < len(dico):
                if seq_codantes[i][j] == clés[k]:
                    tab[i].append(list(dico.values())[k])
                k = k + 1
                
            j = j + 1
            k = 0
            
        i = i + 1
        j = 0
        
    return tab


def adn_encode_molecule(adn, dico, molecule):
    verif = arn(adn)
    if verif == None : 
        return False
    triplés = arn_to_codons(verif)
    
    stop = codons_stop(dico)
    i = 0
    
    while i < len(triplés):
        j = 0
        k = 0
        while j < len(stop):
            if stop[j] in triplés[i]:
                triplés.remove(stop[j])
            j = j + 1
        
        codons = codons_to_aa(triplés, dico)
        
        while k < len(molecule):
            if molecule[k] in codons[i]:
                return True
            k = k + 1
        i = i + 1
        
    return False
