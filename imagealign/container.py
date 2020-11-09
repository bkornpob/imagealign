# Kornpob Bhirombhakdi
# kbhirombhakdi@stsci.edu

import os,glob,warnings

class Container:
    """
    Container is a class facilitating how the output files would be named, and packaged in a folder.
    Name convention is ./savefolder/saveprefix_suffix.ext.
    suffix and ext would be chosen internally in the pipeline. If it is a graphic file such as plots, plotformat determines its ext.
    overwrite, if set True, an existing folder with the same name would be removed at the beginning.
    """
    def __init__(self,saveprefix,savefolder,plotformat,overwrite):
        self.data = {'saveprefix':saveprefix,
                     'savefolder':savefolder,
                     'plotformat':plotformat,
                     'overwrite':overwrite
                    }
        self._prep_container()
    def _prep_container(self):
        tmp = glob.glob('*')
        if self.data['savefolder'] in tmp:
            sentinel = True
        else:
            sentinel = False
        if not sentinel:
            os.mkdir(self.data['savefolder'])
        else:
            if not self.data['overwrite']:
                string = 'Folder {0} already exists. To create a fresh folder, set overwrite = True'.format(self.data['savefolder'])
                warnings.warn(string)
            else:
                os.system('rm -r ./{0}'.format(self.data['savefolder']))
                os.mkdir(self.data['savefolder'])
