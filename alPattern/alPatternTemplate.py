import pymel.core as pm
from alShaders import alShadersTemplate

class AEalPatternTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('space')
        self.addControl('axis')
        self.addControl('shape')
        self.addControl('frequency')
        self.addControl('offset')

        self.addRemapControls()

        self.addControl('color1')
        self.addControl('color2')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()