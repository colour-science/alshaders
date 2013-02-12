from mtoa.ui.ae.shaderTemplate import ShaderAETemplate

class alShadersTemplate(ShaderAETemplate):

    def addRemapControls(self):
        self.beginLayout('Remap', collapse=True)
        self.addControl('RMPinputMin', label='Input min')
        self.addControl('RMPinputMax', label='Input max')

        self.beginLayout('Contrast', collapse=False)
        self.addControl('RMPcontrast', label='Contrast')
        self.addControl('RMPcontrastPivot', label='Pivot')
        self.endLayout() # end Contrast

        self.beginLayout('Bias and gain', collapse=False)
        self.addControl('RMPbias', label='Bias')
        self.addControl('RMPgain', label='Gain')
        self.endLayout() # end Bias and gain

        self.addControl('RMPoutputMin', label='Output min')
        self.addControl('RMPoutputMax', label='Output max')

        self.beginLayout('Clamp', collapse=False)
        self.addControl('RMPclampEnable', label='Enable')
        self.addControl('RMPthreshold', label='Expand')
        self.addControl('RMPclampMin', label='Min')
        self.addControl('RMPclampMax', label='Max')
        self.endLayout() # end Clamp

        self.endLayout() # end Remap

    
    