import pymel.core as pm
from alShaders import alShadersTemplate

class AEalSurfaceTemplate(alShadersTemplate):
    def setup(self):
        self.addSwatch()
        self.beginScrollLayout()

        self.beginLayout('Diffuse', collapse=False)

        self.addControl('diffuseStrength', label='Strength')
        self.addControl('diffuseColor', label='Color')
        self.addControl('diffuseRoughness', label='Roughness')

        self.beginLayout('Subsurface scattering')
        self.addControl('sssMix', label='Mix')
        self.addControl('sssRadius', label='Distance')
        self.addControl('sssRadiusColor', label='Color')
        self.addControl('sssDensityScale', label='Density scale')
        self.endLayout() # end Subsurface scattering

        self.beginLayout('Advanced')
        self.addControl('diffuseExtraSamples', label='Extra samples')
        self.addControl('diffuseEnableCaustics', label='Enable caustics')
        self.endLayout() # end Advanced

        self.endLayout() # end Diffuse

        self.beginLayout('Specular 1', collapse=False)
        self.addControl('specular1Strength', label='Strength')
        self.addControl('specular1Color', label='Color')
        self.addControl('specular1Roughness', label='Roughness')
        self.addControl('specular1Ior', label='IOR')

        self.beginLayout('Advanced')
        self.addControl('specular1RoughnessDepthScale', label='Roughness depth scale')
        self.addControl('specular1ExtraSamples', label='Extra samples')
        self.addControl('specular1Normal', label='Normal')
        self.endLayout() # end Advanced

        self.endLayout() # end Specular 1

        self.beginLayout('Specular 2')
        self.addControl('specular2Strength', label='Strength')
        self.addControl('specular2Color', label='Color')
        self.addControl('specular2Roughness', label='Roughness')
        self.addControl('specular2Ior', label='IOR')

        self.beginLayout('Advanced')
        self.addControl('specular2RoughnessDepthScale', label='Roughness depth scale')
        self.addControl('specular2ExtraSamples', label='Extra samples')
        self.addControl('specular2Normal', label='Normal')
        self.endLayout() # end Advanced
        
        self.endLayout() # end Specular 1

        self.beginLayout('Transmission')
        self.addControl('transmissionStrength', label='Strength')
        self.addControl('transmissionColor', label='Color')
        self.addControl('transmissionLinkToSpecular1', label='Link to spec')
        self.addControl('transmissionRoughness', label='Roughness')
        self.addControl('transmissionIor', label='IOR')

        self.beginLayout('Scattering')
        self.addControl('ssStrength', label='Strength')
        self.addControl('ssTargetColor', label='Color')
        self.addControl('ssDirection', 'Direction')
        self.addControl('ssBalance', label='Balance')
        self.addControl('ssInScattering', label='In-scattering')
        self.addControl('ssDensityScale', 'Density scale')
        self.addControl('ssSpecifyCoefficients', label='Specify coefficients')
        self.addControl('ssScattering', label='Scattering')
        self.addControl('ssAbsorption', label='Absorption')
        self.endLayout() # end Scattering

        self.beginLayout('Advanced', collapse=True)
        self.addControl('transmissionRoughnessDepthScale', label='Roughness depth scale')
        self.addControl('transmissionExtraSamples', label='Extra samples')
        self.addControl('transmissionEnableCaustics', label='Enable internal reflections')
        self.endLayout() # end Advanced

        self.endLayout() # end Transmission

        self.beginLayout('Emission', collapse=True)
        self.addControl('emissionStrength', label='Strength')
        self.addControl('emissionColor', label='Color')
        self.endLayout() # end Emission

        self.beginLayout('IDs', collapse=True)
        for i in range(0,8):
            self.addControl('id'+str(i))
        self.endLayout() # end IDs

        self.beginLayout('Advanced', collapse=True)
        self.addControl('lightGroupsIndirect', label='Indirect light groups')
        self.endLayout() # end Advanced

        self.addBumpLayout()

        pm.mel.AEdependNodeTemplate(self.nodeName)

        self.addExtraControls()
        self.endScrollLayout()