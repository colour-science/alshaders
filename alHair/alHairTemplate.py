import pymel.core as pm
from alShaders import alShadersTemplate

class AEalHairTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.beginLayout('Fibre properties', collapse=False)
        self.addControl('hairColor', label='Color')
        self.addControl('specularWidth', label='Highlight width')
        self.addControl('specularShift', label='Highlight shift')
        self.addControl('opacity')
        self.endLayout() # end Fibre properties

        self.beginLayout('Diffuse', collapse=False)
        self.addControl('diffuseStrength', label='Strength')
        self.addControl('diffuseColor', label='Tint')
        self.addControl('Forward scattering')
        self.addControl('Back scattering')
        self.endLayout() # end Dual scattering

        self.beginLayout('Specular 1', collapse=False)
        self.addControl('specular1Strength', label='Strength')
        self.addControl('specular1Color', label='Tint')
        self.endLayout() # end Specular 1

        self.beginLayout('Specular 2', collapse=False)
        self.addControl('specular2Strength', label='Strength')
        self.addControl('specular2Color', label='Tint')
        self.addControl('glintStrength', label='Glint strength')
        self.addControl('twist', label='Twist')
        self.endLayout() # end Specular 2

        self.beginLayout('Transmission', collapse=False)
        self.addControl('transmissionStrength', label='Strength')
        self.addControl('transmissionColor', label='Tint')
        self.endLayout() # end Transmission

        self.beginLayout('Advanced', collapse=True)
        self.addControl('extraSamples', label='Extra samples')
        self.addControl('transmissionRolloff', label='Transmission rolloff')
        self.addControl('glintRolloff', label='Glint rolloff')
        self.addControl('singleSaturation')
        self.addControl('multipleSaturation')
        self.addControl('dualDepth')
        self.endLayout() # end Advacned

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()