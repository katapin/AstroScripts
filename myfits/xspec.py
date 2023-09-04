"""Extend functionality of pyxspec."""

import os
from dataclasses import dataclass
from contextlib import contextmanager
# from mypython import DEB


class Xspec():
    """Extension for pyxspec."""
    
    __instance = None
    def __new__(cls, *args, **kwargs):
        """Make the class singleton."""
        if cls.__instance:
            return cls.__instance
        else:
            cls.__instance = super().__new__(cls)
            return cls.__instance
    
    def __init__(self, set_rebin=(5,5), statMethod='cstat', show_background=False, xaxis='keV', 
                 query='yes', nTries=5):
        import xspec
        self.__dict__.update(xspec.__dict__)
        with self.__silence():
            self.Plot.xAxis = xaxis
            self.Plot.background=show_background
            self.Plot.setRebin(*set_rebin)
            self.Fit.query=query
            self.Fit.statMethod = statMethod
            self.nTries=nTries
  
    def load_spectra(self, spectra: list , clear=True, quiet=False):
        """Load multiple spectra."""
        #TODO: relative paths of backfile, respfile etc. may cause conflicts with change_dir()
        startWD=os.getcwd()
        startChatter=self.Xset.chatter
        
        if clear:
          self.AllData.clear()
          
        res=[]  #Return list of xspec.Spectrum() objects
        for entry in spectra:
            os.chdir(os.path.dirname(entry.specfile))  #Go the spectrum host dir to load correct bkg and resp files from its FITS-header
            with self.__silence():  #Hide the output because it shows message with no bkg and resp files which may be confusing
                sp = self.Spectrum(entry.specfile)
                if entry.bkgfile:
                    sp.background=entry.bkgfile
                if entry.respfile:
                    sp.response=entry.respfile
                if entry.arffile:
                    sp.response.arf=entry.arffile
                if entry.ignore:
                    for expression in entry.ignore:
                        sp.ignore(expression)
                
                os.chdir(startWD)
            res.append(sp)
            if not quiet:
                sp.show()
    
        return res
    
    def load_single_spectrum(self, specfile:str, bkgfile:str='', respfile:str='', arffile:str='', ignore:list=[], clear=True, quiet=False):
        """Append new spectrum to already loaded."""
        spec = SpecDescription(specfile, bkgfile, respfile, arffile, ignore)
        return self.load_spectra([spec], clear, quiet)[0]
  
    def standard_analysis(self, quiet=False):
        """Fit the model and return parameter values with errors for single spectrum."""
        free_parameters=[]
        for i in range(1,6):
            try:
                mo=self.AllModels(i)
                for k in range(1,mo.nParameters+1):
                    par=mo(k)
                    if not par.frozen: free_parameters.append(k)
            except:
                break
        
        curDelta=100
        curTry=0
        with self.__silence(quiet):
            while (curDelta >= self.Fit.criticalDelta) and (curTry < self.nTries):
                self.Fit.perform()
                first_fit = self.Fit.statistic
                self.Fit.error('1.0 {}'.format(' '.join(map(str,free_parameters))))
                self.Fit.perform()
                second_fit = self.Fit.statistic
                curDelta = first_fit - second_fit
                curTry+=1
            else:
                if curTry == self.nTries:
                    print("Can't fit the model. Tries are exhausted")
                    return None
  
  
    def model_backup(self, index:int=1):
        """Get ModelDescription object."""
        return ModelDescription(self.AllModels(index))
    
    def model_restore(self, model_description):
        """Load model from the ModelDescription object."""
        with self.__silence(True):
            mdl=self.Model(model_description.expression)
            for i in range(model_description.nParameters):
                mdl(i+1).values = model_description[i].values

    def model_absorbtion_remove(self, comp_name:str='TBabs', index:int=1):
        """Remove interstellar absorbtion model component."""
        mds = self.model_backup(index)  #Model description obejct
        self.old_model = mds
        encounters = sum(1 for comp in mds.components if comp.realname == comp_name )
        if encounters == 0:
            raise ValueError(f"Model componen '{comp_name}' is not present in the "\
                    "current model")
        elif encounters > 1:
            raise NotImplementedError("Model component encounters twice or more times. "\
                    "Not supported yet.")
                
        new_expression = mds.expression[len(comp_name):]
        if new_expression[0] == '*': new_expression = new_expression[1:]
        with self.__silence(True):
            mdl=self.Model(new_expression)
            for compname in mdl.componentNames:
                comp = getattr(mdl,compname)       #Model component, e.g. mdl.diskbb
                for parname in comp.parameterNames:
                    par = getattr(comp, parname)   #Parameter  mdl.diskbb.Tin
                    par.values = mds[compname][parname].values
        
    
    @contextmanager
    def __silence(self, quiet=True):
        if quiet:
            startChatter = self.Xset.chatter
            self.Xset.chatter=0
            yield
            self.Xset.chatter=startChatter
        else:
            yield

@dataclass
class SpecDescription():
    specfile:str
    bkgfile:str   = None
    respfile:str  = None
    arffile:str   = None
    ignore:tuple  = None

        
class ModelDescription():
    def __init__(self, model_object):
        if model_object:
            self.expression=model_object.expression
            self.nParameters=model_object.nParameters
            self.__parameters=[]    #Index is xspec's parameter index (from 1)
            self.__components={}    #Key is comp name
            for ind, comp_name in enumerate(model_object.componentNames):
                comp = self.__Component(model_object.__getattribute__(comp_name), ind)
                self.__components[comp_name] = comp
                self.__parameters += comp.parameters
                
    def __getitem__(self, key):
        if type(key) is int:
            return self.__parameters[key]
        elif type(key) is str:
            return self.__components[key]
        else:
            raise TypeError('key must be string or integer')
    
    def __repr__(self):
        return "<ModelDescription '{}' at {}>".format(self.expression, hex(id(self)))
    
    # def __del_comp(self, index):
    #     for par in self.__parameters:
    #         if par.comp_index == index:
    #             self.__parameters.remove(par)
        
    # def del_comp(self, index):
    #     new_instance = copy.deepcopy(self)
    #     new_instance.__del_comp(index)
    #     return new_instance
    
    def print_summary(self):
        print(f"Model '{self.expression}'")
        for par in self.__parameters:
            print(' {:2} {:10}:{:8}\t=\t{}'.format(par.comp.index, par.comp.name, 
                   par.name, str(par.values)))
            
    @property 
    def components(self):
        """Return tuple of model components."""
        return tuple(self.__components.values())
    
    @property 
    def compnames(self):
        """Return tuple of component names."""
        return tuple(self.__components.keys())
    
    @property 
    def parameters(self):
        """Return tuple of model components."""
        return tuple(self.__parameters)
            
    class __Component():
        def __init__(self, comp_object, index: int):
            self.__name=comp_object.name
            self.__index=index
            self.__parameters={}
            for par_name in comp_object.parameterNames:
                self.__parameters[par_name]=self.__Parameter(comp_object.__getattribute__(par_name), self)
                
        def __getitem__(self, key):
            if not type(key) is str:
                raise TypeError('index must be string')
            return self.__parameters[key]
        
        def __repr__(self):
            return "<Component '{}' at {}>".format(self.name, hex(id(self)))
        
        @property
        def name(self):
            return self.__name
        
        @property
        def realname(self):
            return self.__name.split('_')[0]
        
        @property
        def index(self):
            return self.__index
            
        @property
        def parameters(self):
            return list(self.__parameters.values())

        @property
        def parameter_names(self):
            return list(self.__parameters.keys())

        
        class __Parameter():
            def __init__(self, param_object, comp):
                self.__name=param_object.name
                self.__index=param_object.index   #Xspec index, starting from 1
                self.__comp = comp
                self.values=param_object.values
                
            def __repr__(self):
                return "<Parameter {}={} of '{}' at {}>".format(self.name, self.values[0],
                    self.__comp.name, hex(id(self)))
            
            @property
            def name(self):
                return self.__name
            
            @property 
            def index(self):
                return self.__index
            
            @property 
            def comp(self):
                return self.__comp 
            
            
#### Modify chain ####

class _ChainColInfo():
    """Class to store properties of xspec's chain FTIS tables."""
    
    def __init__(self, colobj):
        """
        Parameters
        ----------
        colobj : astropy.io.fits.column.Column
        """
        pair=colobj.name.split('__')
        self.parname = pair[0]
        self.index = int(pair[1])   #Real pararmter's index, not column number
        self.longname = colobj.name
        self.unit = colobj.unit

def modify_chain(infile, outfile, params_to_remove:list):
    """Remove columns from the chain filo mimic different model.
    
    Xspec's chain files contain columns with naming 'parname__index',
    (for example, PhoIndex__2). It's important to note that the index here 
    refers not to indexes of columns in the chain FITS table but to
    native xspec paramter indixes. Since frozen parameters are absent 
    in the chain tables, this indexind is uneven. There for 'params_to_remove'
    must list indixes of all the paramters of the model componen beeing
    removed, included the fozen ones (and thus absent in the chain table).
    """
    from astropy.io import fits
    from astropy.time import Time
    from datetime import datetime
    
    #TODO: create the params_to_remove list from xspec.Model
    if not isinstance(params_to_remove,list):
        raise TypeError("Wrong type of the 'param_to_remove' argument")
  
    with fits.open(infile) as fts:
        tbl=fts[1]
        colinfo = [_ChainColInfo(col) for col in tbl.columns if col.name != 'FIT_STATISTIC']
                
        minpar = min(params_to_remove)  #Minimum index of parameters to remove
        for col in colinfo:
            if col.index < minpar:
                continue
            elif col.index in params_to_remove:
                tbl.columns.del_col(col.longname)
            else:
                #Count removed parameters staying before the current column
                removed_before = sum(1 for x in params_to_remove if x < col.index)
                newindex = col.index - removed_before
                newname = col.parname+'__'+str(newindex)
                tbl.columns.change_name(col.longname, newname) 
                
                #Bugfix: edit wrong units of renamed columns
                #Otherwise xspec rejects this file and fv crashes
                tbl.columns.change_attrib(newname, 'unit', col.unit) 
            
        now = Time(datetime.now())
        tbl.header.update()
        
        #Bugfix: remove TUNITn if 'unit' is None
        for i in range(1, len(tbl.columns)+1): 
            curTUNIT = f'TUNIT{i}'
            if curTUNIT in tbl.header:
                if tbl.header[curTUNIT] is None:
                    del tbl.header[curTUNIT]
        tbl.header.add_history(f'Produced from {infile} on {now.fits}')
        print(tbl)
        tbl.writeto(outfile)