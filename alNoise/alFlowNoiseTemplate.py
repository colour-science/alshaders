import pymel.core as pm
from alShaders import alShadersTemplate

class AEalFlowNoiseTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('space')
        self.addControl('frequency')
        self.addControl('octaves')
        self.addControl('lacunarity')
        self.addControl('gain')
        self.addControl('angle')
        self.addControl('advection')
        self.addControl('turbulent')

        self.addRemapControls()

        self.addControl('color1')
        self.addControl('color2')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()