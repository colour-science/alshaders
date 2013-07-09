def makeAlSurfaceAOVs():
	standardCompatibleAOVNames = [
		"Arnold_Direct_Diffuse",
		"Arnold_Indirect_Diffuse",
		"Arnold_Direct_Specular",
		"Arnold_Indirect_Specular",
		"Arnold_SSS",
		"Arnold_Refraction",
		"Arnold_Emission",
		]
		
	alshadersAOVs = [
		"Arnold_Diffuse_Color",
		"Arnold_Direct_Backlight",
		"Arnold_Direct_Diffuse_raw",
		"Arnold_Indirect_Backlight",
		"Arnold_Indirect_Diffuse_raw",
		"Arnold_Direct_Specular_2",
		"Arnold_Indirect_Specular_2",
		"Arnold_Single_Scatter",
		"Arnold_uv",
		"Arnold_depth",
		]
		
	alshadersLightgroupsAOVs = [
		"light_group_1",
		"light_group_2",
		"light_group_3",
		"light_group_4",
		"light_group_5",
		"light_group_6",
		"light_group_7",
		"light_group_8",
		]
		
	alShadersIDAOVs = [
		"id_1",
		"id_2",
		"id_3",
		"id_4",
		"id_5",
		"id_6",
		"id_7",
		"id_8",
		]
		
	aovList = standardCompatibleAOVNames + alshadersAOVs + alshadersLightgroupsAOVs + alShadersIDAOVs

	channelsContainer = Application.Dictionary.GetObject("Passes.RenderOptions.Channels", False)
	channelList = []

	for channel in channelsContainer.NestedObjects:
		channelList.append(channel.name)

	for aov in aovList:
		if aov not in channelList:
			Application.LogMessage("Adding: %s" % aov)
			Application.CreateRenderChannel(aov, "siRenderChannelColorType", "")

makeAlSurfaceAOVs()