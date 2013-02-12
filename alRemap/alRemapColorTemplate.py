import pymel.core as pm
from alShaders import alShadersTemplate

class AEalRemapColorTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('input')
        self.addControl('gamma')
        self.addControl('saturation')
        self.addControl('hueOffset')

        self.beginLayout('Contrast', collapse=False)
        self.addControl('contrast')
        self.addControl('contrastPivot', label='Pivot')
        self.endLayout() # end Contrast

        self.addControl('gain')
        self.addControl('exposure')
        self.addControl('mask')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()