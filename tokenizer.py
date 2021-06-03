from random import shuffle
from rdkit import Chem
from molvs.standardize import Standardizer



class ChemTokenizer:

    def __init__(self, atoms_list: list = None, other_atoms_list: list = None):

        atoms = [{'Ac', 'Ag', 'Al', 'Am', 'Ar', 'As', 'At', 'Au', 'Ba', 'Be', 'Bi', 'Bk', 'Br', 'Ca', 'Cd',
                  'Ce', 'Cf', 'Cl', 'Cm', 'Cr', 'Cu', 'Dy', 'Er', 'Es', 'Eu', 'Fe', 'Fm', 'Fr', 'Ga', 'Gd',
                  'Ge', 'He', 'Hf', 'Hg', 'Ir', 'Kr', 'La', 'Li', 'Lr', 'Lu', 'Md', 'Mg', 'Mn', 'Mo', 'Na',
                  'Nd', 'Ne', 'Ni', 'Pa', 'Pd', 'Pm', 'Pr', 'Pt', 'Pu', 'Ra', 'Rb', 'Re', 'Rh', 'Ru', 'Se',
                  'Si', 'Sm', 'Sr', 'Ta', 'Tb', 'Tc', 'Te', 'Th', 'Ti', 'Tl', 'Tm', 'Xe', 'Zn', 'Zr', 'in',
                  'se', 'te'},
                 
                 {'#', '$', '%', '(', ')', '*', '+', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6',
                  '7', '8', '9', ':', '=', '@', 'B', 'C', 'D', 'F', 'H', 'I', 'K', 'N', 'O', 'P', 'S',
                  'U', 'V', 'W', 'Y', '[', '\\', ']', 'b', 'c', 'n', 'o', 'p', 's'}]

        complicated = {'Co', 'Cs', 'Ho', 'In', 'Nb', 'No', 'Np', 'Os', 'Pb', 'Po',
                       'Sb', 'Sc', 'Sn', 'Yb', 'os', 'pb', 'sb', 'sc', 'sn'}

        if atoms_list:
            self.atoms = atoms_list
        else:
            self.atoms = atoms

        if other_atoms_list:
            self.other = other_atoms_list
        else:
            self.other = complicated

    @staticmethod
    def _to_check(result: 'tokenized SMILES', SMILES: 'src SMILES') -> bool:
        
        """
        Function checks whether tokenization is correct
        
        """
        
        if not ''.join(result) != SMILES:
            return False


    def tokenize(self, SMILES_list: list, check=False) -> list:
        
        """
        Tokenization of SMILES sequences from src list
        
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
    def randomize(SMILES_list: list, random_number: int = 30, max_iter: int = 500) -> list:

        """
        Randomization of each SMILES sequence in the list 
        in order to generate several SMILES representations for the same molecule
        
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
    def standardizer(SMILES_list: list) -> list:

        """
        Standardization of each SMILES in the list using MolVS package
        
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




    def process(self, SMILES_list: list, random_number: int = 30, max_iter: int = 500, check=False) -> list:

        """
        A serial of actions over list of SMILES:
        1) standartization: list of standardized SMILES
        2) randomization: list of lists of RandomSMILES
        3) tokenization: list of lists of tokenized SMILES
        
        """

        SMILES_list = self.standardizer(SMILES_list)
        SMILES_list = self.randomize(SMILES_list, random_number, max_iter)

        new_SMILES = []
        for SMILES in SMILES_list:
            new_SMILES.append(self.tokenize(SMILES, check))

        return new_SMILES
