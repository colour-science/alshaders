import pymel.core as pm
from alShaders import alShadersTemplate

class AEalBlackbodyTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('temperature')
        self.addControl('strength')

        self.beginLayout('Advanced', collapse=True)
        self.addControl('physicalIntensity', label='Physical intensity')
        self.addControl('physicalExposure', label='Physical exposure')
        self.endLayout() # end Advanced

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()