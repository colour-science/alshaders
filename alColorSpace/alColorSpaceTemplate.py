import pymel.core as pm
from alShaders import alShadersTemplate

class AEalColorSpaceTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('input', label='Input')
        self.addControl('sourceSpace', label='Source space')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()