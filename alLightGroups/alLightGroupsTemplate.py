import pymel.core as pm
from alShaders import alShadersTemplate

class AEalLightGroupsTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('color')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()