from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

env.io.atom_files_directory = ['.', '../atom_files']

a = automodel(env, alnfile = 'alignment.ali',knowns = 'mT1', sequence = 'hT1')

a.starting_model= 1
a.ending_model  = 1

a.make()
