import pymel.core as pm
from alShaders import alShadersTemplate

class AEalFractalTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.addControl('mode')
        self.addControl('space')
        self.addControl('scale')
        self.addControl('frequency')
        self.addControl('time')
        self.addControl('octaves')
        self.addControl('distortion')
        self.addControl('lacunarity')
        self.addControl('gain')
        self.addControl('turbulent')
        self.addControl('ridged')
        self.addControl('ridgeOffset')

        self.addRemapControls()

        self.addControl('color1')
        self.addControl('color2')

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()
