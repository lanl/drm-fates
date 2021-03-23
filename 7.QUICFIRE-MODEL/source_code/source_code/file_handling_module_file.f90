	MODULE file_handling_module
	
	implicit none
	
	character(256) :: workingDirectory
	character :: fileSeparator
	
	integer,parameter ::												&
		STR_LEN = 500,													&
		!STDOUT 											= 6,			&
		!DBGFILE											= 10,			&
		ID_FILE_CPU										= 50,			&
		ID_FILE_Z_QU									= 60,			&
		ID_FILE_BLDG_ADVANCED						= 80,			&
		ID_FILE_QFIRE_INP								= 81,			&
		ID_FILE_GRIDLIST								= 82,			&
		ID_FILE_QF_FUELDENSITY						= 83,			&
		ID_FILE_TREESRHOF								= 84,			&
		ID_FILE_QF_FUELMOIST							= 85,			&
		ID_FILE_TREESMOIST							= 86,			&
		ID_FILE_QF_FUELHEIGHT						= 87,			&
		ID_FILE_TREESFUELDEPTH						= 88,			&
		ID_FILE_QF_IGNITION							= 89,			&		
		ID_FILE_TIMELOG								= 91,			&
		ID_FILE_TRAJ									= 101,		&
		ID_FILE_TRAJ_MERGE							= 102,		&
		ID_FILE_WIND_QU_U								= 110,		&
		ID_FILE_WIND_QU_V								= 111,		&
		ID_FILE_WIND_QU_W								= 112,		&		
		ID_FILE_PM_EMISSIONS_DISTR					= 130,		&
		ID_FILE_PM_EMISSIONS							= 131,		&
		ID_FILE_CO_EMISSIONS							= 132,		&
		ID_FILE_FUEL_OUT								= 133,		&
		ID_FILE_WIND_FIRE_U							= 134,		&
		ID_FILE_WIND_FIRE_V							= 135,		&
		ID_FILE_WIND_FIRE_W							= 136,		&
		ID_FILE_WIND_FIRE_SIGMA						= 137,		&
		ID_FILE_MOISTURE_OUT							= 139,		&
		ID_FILE_FIREINDEX								= 140,		&
		ID_FILE_FUELHEIGHT							= 141,		&
		ID_FILE_MASSBURNT								= 144,		&
		ID_FILE_REACTRATE								= 145,		&
		ID_FILE_EN2ATM									= 146,		&
		ID_FILE_CONVHUMAN								= 147,		&
		ID_FILE_THERMDOSE								= 148,		&
		ID_FILE_IGNPATTERN							= 150,		&
		ID_FILE_RASTERORIGIN							= 152,		&
		ID_FILE_IGNITE_SEL							= 153,		&
		ID_FILE_IGNITE									= 154,		&
		ID_FILE_FB_OUT									= 160,		&
		ID_FILE_FB_IN									= 161,		&
		ID_FILE_FB_LOST								= 162,		&
		ID_FILE_ADV_PLUME								= 170,		&
		ID_FILE_QU_SIMPARAMS							= 207,		&
		ID_FILE_QU_BLDFLAG							= 208,		&
		ID_FILE_QU_BLDARRAY							= 209,		&		
		ID_FILE_QP_SOURCE								= 217,		&		
		ID_FILE_QU_BLD									= 225,		&
		ID_FILE_QU_LANDUSE							= 226,		&		
		ID_FILE_QU_METPARAMS							= 228,		&
		ID_FILE_QU_FILEOPT							= 229,		&		
		ID_FILE_QU_CELLTYPE_BIN						= 250,		&
		ID_FILE_QU_BUILDOUT							= 251,		&
		ID_FILE_QU_VELOCITY_BIN						= 252,		&
		ID_FILE_QU_UNDIST_BIN						= 253,		&
		ID_FILE_QU_VEG_BIN							= 255,		&		
		ID_FILE_QU_CELLTYPE_DAT						= 260,		&
		ID_FILE_QU_VELOCITY							= 262,		&
		ID_FILE_QU_SURFACE							= 318,		&
		ID_FILE_SENSOR									= 319,		&
		ID_FILE_QU_UOFIELD							= 320,		&
		ID_FILE_QU_STAGGERED							= 321,		&
		ID_FILE_QU_PARAMRANGE						= 322,		&
		ID_FILE_QU_ERROR								= 325,		&
		ID_FILE_QU_MOVINGCOORD						= 363
	
	!integer :: msgoutfile

	END MODULE file_handling_module
!===========================================================================	
!===========================================================================
