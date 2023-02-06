"""
Response Function Library
Primary Author: Paul O'Brien, The Aerospace Corporation paul.obrien@aero.org
"""

class ChannelResponse(object):
    """
    Virtual base class for channel responses
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,RESP_TYPE=None,**kwargs):
        """
        channel = ChannelResponse.create(resp_type,...)
            RESP_TYPE is reponse type '[E]',<etc>,'[E,TH,PH]'
        channel = ChannelResponse.create(info,...)
            info is a dict of response type data (e.g., loaded from a file)
        returns an initialized response object of the specified type
        """
        if isinstance(RESP_TYPE,cls): # undocumented feature
            raise Warning("Called ChannelReponse.create on object. Use object's copy method instead")
            return cls.copy()
        if isinstance(RESP_TYPE,dict):
            kwargs = {**RESP_TYPE,**kwargs} # merge response type and kwargs
            RESP_TYPE = kwargs['RESP_TYPE']
        assert isinstance(RESP_TYPE,str),'response type (RESP_TYPE) can only be string'
        RESP_TYPE = RESP_TYPE.upper()
        if RESP_TYPE == '[E]':
            return CR_Omni.create(**kwargs)
        elif RESP_TYPE == '[E][TH]':
            return CR_Esep_sym.create(**kwargs)
        elif RESP_TYPE == '[E][TH,PHI]':
            return CR_Esep_Asym.create(**kwargs)
        elif RESP_TYPE == '[E,TH]':
            return CR_insep_sym.create(**kwargs)
        elif RESP_TYPE == '[E,TH,PH]':
            return CR_insep_Asym.create(**kwargs)
        else:
            raise Exception('Unknown response type (RESP_TYPE): %s' % RESP_TYPE)
    def copy(self):
        """
        channel.copy() return a copy of this object
        """
        ChannelResponse.create(**self.to_dict())
    def to_dict(self):
        """
        channel.to_dict() return a nested dict of this object
        """
        raise NotImplementedError

class CR_Omni(ChannelResponse):
    """
    CR_Omni
    omnidirectional channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        return cls(**kwargs)
    
class CR_Esep_sym(CR_Omni):
    """
    CR_Esep_sym
    energy-separable, phi-symmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        return cls(**kwargs)    
    
class CR_Esep_Asym(CR_Esep_sym):
    """
    CR_Esep_Asym
    energy-separable, phi-asymmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        return cls(**kwargs)    
    
class CR_insep_sym(CR_Omni):
    """
    CR_insep_sym
    energy-inseparable, phi-symmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        return cls(**kwargs)    
    
class CR_insep_Asym(CR_Omni):
    """
    CR_insep_Asym
    energy-inseparable, phi-asymmetric channel response
    Use factory class method .create() to create
    appropriate subclass from response dictionary
    """
    @classmethod
    def create(cls,**kwargs):
        return cls(**kwargs)    

if __name__ == '__main__':
    # very simple test
    print(ChannelResponse.create(RESP_TYPE='[E,TH]'))
