import pymel.core as pm
from alShaders import alShadersTemplate

class AEalRemapFloatTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('input')

        self.addRemapControls()

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()