import pymel.core as pm
from alShaders import alShadersTemplate

class AEalCurvatureTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('samples', label='Samples')
        self.addControl('sampleOffset', label='Offset')
        self.addControl('sampleRadius', label='Radius')

        self.addRemapControls()

        self.addControl('color1', label='Color 1')
        self.addControl('color2', label='Color 2')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()