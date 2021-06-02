from random import shuffle
from rdkit import Chem
from molvs.standardize import Standardizer



class ChemTokenizer:

    def __init__(self, atoms_list=None, other_atoms_list=None):

        atoms = [{'Ac', 'Bi', 'Be', 'Ta', 'Ba', 'Nd', 'Gd', 'Mg', 'Au', 'Th', 'Na', 'Cr', 'Hf',
                 'Te', 'Mo', 'Pd', 'Eu', 'Pu', 'Cl', 'Ne', 'Mn', 'Zn', 'Sm', 'Hg', 'Re', 'se',
                  'te', 'As', 'Ge', 'Ce', 'Ag', 'Tb', 'Zr', 'Ti', 'Tl', 'Xe', 'Si', 'Cu', 'Sr',
                  'Er', 'He', 'Pr', 'Pt', 'Ga', 'Cd', 'Fe', 'Rh', 'Br', 'At', 'Al', 'Li', 'Ir',
                  'Cm', 'Dy', 'Ca', 'Am', 'Ru', 'Pa', 'Tc', 'Se', 'Tm', 'Lu', 'La', 'Ni', 'in',
                  'Ar', 'Kr'},
                 {'C', '-', 'V', 'p', 'I', 'B', '1', 'U', '8', 'N', 'b', 'P', 's', 'O', '4', 'Y', '#',
                  'n', 'H', 'D', '/', '\\', 'o', 'K', '7', 'F', 'c', '[', ']', '+', '.', ')', '0', '5',
                  '9', 'S', '2', '(', '6', '%', '=', '3', 'W', '@', '*', '$', ':'}]

        complicated = {'Yb', 'Ho', 'Po', 'Sb', 'sb', 'Rb', 'Np', 'Nb', 'Sc', 'sc', 'Pb',
                       'pb', 'Sn', 'sn', 'Co', 'Cs', 'In', 'Os', 'os', 'No'}

        if atoms_list:
            self.atoms = atoms_list
        else:
            self.atoms = atoms

        if other_atoms_list:
            self.other = other_atoms_list
        else:
            self.other = complicated

    @staticmethod
    def _to_check(result, SMILES):
        
        """
        :param: tokenized SMILES
        :param: src SMILES sequence
        
        """
        
        if not ''.join(result) != SMILES:
            return False


    def tokenize(self, SMILES_list, check = False):
        
        """
        :param: src list of SMILES sequences
        :param: whether to check tokenized SMILES or not
        
        """

        new_SMILES = []
        for SMILES in SMILES_list:
            
            result = []
            j = 0
            while j <= len(SMILES)-1:
                if SMILES[j] == '[':
                    if SMILES[j+1:j+3] in self.other:
                        result.append('[')
                        result.append(SMILES[j+1:j+3])
                        j+=3
                    else:
                        result.append('[')
                        j+=1

                elif SMILES[j:j+2] in self.atoms[0]:
                    result.append(SMILES[j:j+2])
                    j+=2
                elif SMILES[j] in self.atoms[1]:
                    result.append(SMILES[j])
                    j+=1
                else:
                    print(f"There's problem with SMILES: {SMILES}")
                    new_SMILES.append(None)
                    break

            if check:
                if self._to_check(result, SMILES):
                    print(f"There's problem with SMILES: {SMILES}")

            new_SMILES.append(result)
                
        return new_SMILES
    

    @staticmethod
    def randomize(SMILES_list, random_number = 30, max_iter = 500):

        """
        :param: src list of SMILES sequences
        :param: number how many random SMILES will be generated
        :param: number how many iterations will be
        
        """

        new_SMILES = []
        for SMILES in SMILES_list:
            
            mol = Chem.MolFromSmiles(SMILES)
            if not mol:
                print(f"SMILES {SMILES} can't be converted into mol format!")
                new_SMILES.append(None)
                continue

            atoms_number = list(range(mol.GetNumAtoms()))
            random_smiles = []
            i = 0
            while len(random_smiles) < random_number and i < max_iter:
                    
                shuffle(atoms_number)
                new_atoms_number = Chem.RenumberAtoms(mol, atoms_number)
                new_smiles = Chem.MolToSmiles(new_atoms_number, canonical = False)
                    
                if new_smiles not in random_smiles:
                    random_smiles.append(new_smiles)
                i += 1
            new_SMILES.append(random_smiles)

        return new_SMILES



    @staticmethod
    def standardizer(SMILES_list):

        """
        :param: src list of SMILES sequences
        :param: number how many random SMILES will be generated
        :param: number how many iterations will be
        
        """

        standard = Standardizer()

        new_SMILES = []
        for SMILES in SMILES_list:

            mol = standard.standardize(Chem.MolFromSmiles(SMILES))
            mol = standard.charge_parent(mol, skip_standardize=True)
            mol = standard.stereo_parent(mol, skip_standardize=True)
            mol = standard.isotope_parent(mol, skip_standardize=True)
            mol = standard.tautomer_parent(mol, skip_standardize=True)

            new_SMILES.append(Chem.MolToSmiles(standard.standardize(mol)))
                
        return new_SMILES




    def process(self, SMILES_list, random_number = 30, max_iter = 500, check = False):

        """
        :param: src list of SMILES sequences
        
        """

        SMILES_list = self.standardizer(SMILES_list)
        SMILES_list = self.randomize(SMILES_list, random_number, max_iter)

        new_SMILES = []
        for SMILES in SMILES_list:
            new_SMILES.append(self.tokenize(SMILES, check))

        return new_SMILES
