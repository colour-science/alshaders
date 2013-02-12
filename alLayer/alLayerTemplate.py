import pymel.core as pm
from alShaders import alShadersTemplate

class AEalLayerTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('layer1', label='Layer 1')
        self.addControl('layer2', label='Layer 2')
        self.addControl('mix', label='Mix')
        self.addControl('debug', label='Debug')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()