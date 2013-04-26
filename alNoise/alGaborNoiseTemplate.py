import pymel.core as pm
from alShaders import alShadersTemplate

class AEalGaborNoiseTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('space')
        self.addControl('frequency')
        self.addControl('anisotropy')
        self.addControl('anistropyDirection')
        self.addControl('filter')
        self.addControl('bandwidth')
        self.addControl('impulses')


        self.addRemapControls()

        self.addControl('color1')
        self.addControl('color2')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()