% This is a configuration setting files that shared among all matlab files.
% This file is to facilitate all the environment-specific variables that
% could be used so that all .m files have no hardcoded path.

%% Pre HCTSA
BLOCKEDFLOAD_DIR='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/library/blockedfloader';
DATA_DIR = '/Volumes/Untitled/edfs';
EDF_FILE = 'ccshs-trec-1800001.edf';
WHICH_DATA = 1;

%% Shared (Common)
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
NUM_CHANNELS = 3;

%% Post HCTSA
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
HCTSA_DATA_DIR='./180001_HCTSA';

%% selectdata and crossval
ANSWER_FILE='ccshs_1800001_annot.mat';
CONFUSION_MATRIX_FILE='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/180001_HCTSA';
CM_SAVE_DIR='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/180001_HCTSA';
HCTSA_FILE='180001_HCTSA/HCTSA_N.mat';
OPS_FILE='reduced_ops.txt';
TRAINING_PERCENTAGE = 0.7;
NUM_CHANNELS_USED_FOR_CROSSVAL = 1;
PLOT_CONFUSION_MATRIX=true;
PLOT_ACCURACY_REPORT=true;

%% Experiment configurations to run
EXPS_TO_RUN=[1:10];
THREAD_POOL=4;
EVAL_THRESHOLD=0.05;
FILTERED_FEATURE_IDS=[5878,5872,5874,5435,2329,2309,5868,5870,872,389,5879,5877,5446,873,3460,874,5436,5372,1859,4047,1825,5880,4262,3380,1853,5734,5882,5447,1880,5493,3431,3405,1849,3070,1661,3355,4302,1854,5450,5439,3390,3375,307,4316,5378,3777,2497,3430,5416,3435,3450,371,5415,3409,2482,370,3136,3139,5955,3388,392,3410,3357,5196,3059,1013,1850,1875,4048,5178,1580,1735,5274,3376,1739,3452,3360,3102,1874,3151,5373,3277,3148,1927,3900,3373,372,4263,2685,5375,2481,5218,1911,3142,1878,5191,2145,3265,3137,5369,420,3134,2499,3103,2143,3237,870,3154,2684,1943,5221,415,364,5226,3966,5856,1851,2798,5280,3157,3433,3145,3354,1569,423,5786,3492,5371,5195,2674,3378,3483,5735,3437,3449,2142,393,3459,2695,3403,5275,1041,3370,3352,5438,5876,5875,3381,3391,2744,5417,3292,3276,5413,5419,369,1981,3612,3769,5483,3472,3841,2140,3461,2799,3699,5732,3465,5884,5486,354,3695,5381,5449,3407,1685,3623,5155,2683,5273,5881,5145,1467,3158,968,3621,3622,5150,2790,3789,3475,5967,1979,4261,2789,1683,1571,5175,869,1737,3253,2781,3818,2114,5279,5197,1463,5663,3453,3822,5349,1980,3182,5200,1855,1881,3379,3865,969,2800,425,3415,3100,5950,436,3610,3611,3852,1882,3099,5485,3308,2675,1398,1681,3968,5199,2771,3436,5710,5225,5379,5272,3404,5102,2056,5370,424,5220,5689,970,3671,2742,3760,3759,5853,434,971,378,1687,5517,428,5604,2780,2693,5794,972,3211,3385,3312,1026,3214,5618,3141,5277,868,5441,3508,1578,5537,3698,5481,5671,5622,5057,1353,3323,4245,3853,1848,5724,5198,5154,5873,5521,2727,3455,5384,5420,5612,5783,3742,450,390,3892,5177,5541,429,3729,5626,26,5347,5939,3439,5597,973,373,440,2764,5471,3874,3397,5110,4329,3335,199,5350,28,2772,5390,5489,384,3505,30,451,3396,439,3779,3288,5414,1461,3443,3901,5176,5395,2705,5201,5701,5376,3387,5578,4942,5156,1873,3702,3179,5601,5422,5172,5245,5945,3328,3152,427,5323,5607,5276,5252,5704,367,1016,3912,5516,376,5520,3147,1204,5620,1059,5476,3425,5557,3751,981,5215,203,5864,3778,1464,5495,3697,1352,5183,1668,5703,1132,5883,3502,5730,3444,2754,5257,2483,976,3495,2720,5223,5849,974,2174,365,5166,2703,3428,3568,442,984,5804,5181,3485,987,3825,5779,5712,5736,3412,2726,382,975,5540,3883,3882,1354,5657,3133,449,3984,5655,2288,5412,3476,2162,24,2540,5055,990,2652,5561,5556,5536,4260,3821,2774,2596,5690,3479,5580,5650,3579,3718,5859,443,5451,2316,1689,5765,5798,2743,5180,5484,3893,426,5202,5702,3791,5337,3351,5801,3847,5267,5130,5560,432,3451,91,3770,5318,3787,132,2694,5695,4042,5194,5666,5799,3910,2719,5452,5248,5848,5728,4083,1266,3135,5705,3567,3566,5265,3324,5718,3140,448,377,5568,1001,2462,2297,5351,5269,306,3468,2108,431,430,3383,5458,3913,3656,5229,5715,3601,3864,3333,3898,5829,3389,2366,5203,3790,4885,5209,2122,5440,305,1679,5711,3903,3097,2752,3915,3836,3902,2463,5837,1264,5222,5173,2151,5835,3138,5382,2495,3457,309,5838,3911,3275,3395,5352,4049,5466,2721,5334,3484,2464,4084,3869,447,3311,5242,3132,982,2387,1675,1680,1919,5322,5828,3497,3605,33,5478,3666,3665,2346,3401,3746,3627,3487,2465,3118,5394,453,5423,5473,2141,3119,35,5731,5951,5649,1031,4265,1568,5335,3445,2493,5789,2072,3887,3820,5332,965,308,3780,395,3599,3600,3442,3788,452,3827,3856,1879,1810,2080,2222,5717,5425,3096,5824,5576,3200,5598,908,2236,5582,3478,3741,3094,5429,3107,3197,2713,1674,5165,5164,3673,3730,4043,5833,1935,2457,967,964,5426,2229,5572,3792,5532,5213,5333,5533,5726,2775,3578,3577,3667,3583,5496,3272,3616,916,214,1116,2381,3914,1671,4044,366,694,1251,1466,1912,4957,988,1951,4045,909,5162,1647,3947,2326,5564,5500,1579,3713,5760,5227,5161,5160,409,5259,966,3859,5109,3572,2792,4612,3955,3454,2458,4031,5595,2498,5491,3202,3626,5428,5803,5163,3305,3585,3704,2267,2770,2491,3474,2783,5570,3199,2375,5672,5316,3350,5738,202,1344,1823,3625,2492,3342,454,985,3895,1053,3908,1198,381,1249,1567,3764,5847,3668,21,5432,3629,37,1126,5494,2682,963,3574,5836,4011,5437,5769,5709,459,1612,961,5508,210,5219,1667,1672,5147,3890,5249,2376,5179,3842,5317,5217,215,3507,2732,5434,209,1794,2352,4082,5694,5687,230,5588,3315,3341,2477,2490,255,5627,4869,5465,4249,1043,5261,2461,3274,5867,259,250,1263,1248,5159,2227,3470,5260,2459,962,3860,991,1734,363,3504,2064,1404,3736,2170,5851,2389,960,5716,5129,387,3259,20,2264,3839,2276,229,2334,3300,3332,5569,5146,3905,5885,4233,2027,5676,5468,1382,3069,5336,2762,6025,5674,2149,3772,3155,205,379,5675,356,5480,39,1256,1341,5733,1021,926,394,3615,910,4682,3093,5211,5629,2374,3130,5630,5385,5388,3129,4927,5806,598,2339,693,3440,2793,2161,5186,2784,5631,3737,3614,3607,3256,3245,6030,2277,2363,260,313,2296,1413,1188,5947,5793,1362,3849,1167,3153,5477,3967,3325,383,5253,1686,3400,3654,3655,3618,3501,5236,2535,4574,5575,5448,2275,5321,5510,1203,5571,5497,5152,3106,3696,5472,445,4062,2332,4060,865,5563,3782,3843,2343,3675,2394,2287,3917,5456,4021,5573,864,1666,2284,5205,1202,3670,1236,254,2647,2422,5808,4088,5778,3427,1403,5362,2591,5942,456,2173,3979,3669,2395,5149,3128,2396,5937,5467,3833,2371,1622,1942,1630,31,5809,2489,90,5962,5331,5810,1926,5719,1632,1414,5780,5401,3168,1255,5247,5559,5740,2115,3231,3441,2260,4628,5501,5445,3794,813,3959,5742,3165,483,5741,2769,5240,824,820,5206,1738,5869,2207,5700,3464,3372,5153,3170,4310,1665,3851,996,3399,3167,994,2295,4893,225,1241,2741,5207,3358,2335,2460,1617,5865,1343,5767,3733,3326,1676,2110,5479,3348,999,5850,4019,2415,5474,5512,5519,599,5459,5528,4655,5936,3344,5661,5462,866,4254,2430,5613,2315,2292,5410,5230,3221,2228,2157,6029,3960,5357,1259,127,1130,5158,2306,2456,2182,1688,5935,692,2202,2307,2973,5470,863,950,441,5363,2211,5356,5522,4894,3488,5125,3250,5802,4085,5400,3467,2159,3466,5444,92,3204,5857,139,5463,5433,2403,4623,5761,5105,128,5475,2311,5270,3413,5659,408,5455,5934,2303,126,5348,2111,4041,2440,2232,2154,1556,5238,93,5553,1004,5953,5523,2677,374,3848,224,1788,407,2167,5174,4257,5053,840,3884,3886,5460,1673,5513,5046,3281,881,2112,1131,3346,3196,5552,437,3131,3058,3757,1244,1852,2327,2158,3226,3251,2057,129,3628,2242,1909,1239,2673,1242,5871,5461,2230,5324,5545,2383,2078,1057,4678,2331,5502,1006,5843,138,2225,1338,22,125,137,959,3203,4892,3205,2437,1096,2704,5933,438,1910,980,4677,3885,201,2678,204,1925,5457,458,5151,5565,5544,3295,380,5360,5359,375,2063,3318,2951,2712,455,1237,4604,4871,2165,3319,2808,1941,1058,2055,1809,3220,3247,5518,3463,220,3637,5103,5171,5574,1333,3880,3195,5771,3636,5766,5224,3819,1252,3338,2426,249,4055,2488,314,2147,5214,939,3017,5047,2442,5584,1658,200,1635,3236,3233,3838,3482,2848,5393,3763,839,5271,2435,983,5354,4890,3228,3837,2323,5890,3970,4606,4607,5498,989,5320,2450,5144,1231,4033,2447,1232,3386,2272,5682,3434,2336,3055,5983,2289,1702,5148,2224,310,2449,3302,5314,444,3184,2444,1790,190,2455,986,2895,5611,4959,1599,4237,2391,5509,1678,3398,4246,3207,19,4605,4903,706,189,1877,2718,5845,2692,124,3471,5844,1598,2235,1939,3582,2239,1871,705,2302,2428,482,5319,3581,3267,3493,3680,2850,2766,3186,2280,2397,2669,1332,5619,29,2733,5353,1876,5529,4061,3382,3232,1233,5123,5366,5104,2438,5392,3481,2355,1664,1660,2431,5645,1663,1793,2321,2031,3218,385,3216,2714,2109,1659,5535,3369,1165,2401,5469,4609,219,2995,1834,5814,4022,2443,1407,1558,4650,4012,2237,5232,410,4087,3624,5756,5562,94,3175,1166,5800,3377,5841,1260,3039,2496,5499,5623,1696,40,4032,1199,3115,3095,3961,3958,879,2048,1157,5691,5889,264,3240,3406,3719,5157,1153,5713,1152,2234,2433,2359,3209,136,2351,1002,5932,1161,3260,2166,1397,5548,3710,4267,86,27,2345,3368,1590,4603,123,1331,1163,1865,3314,2436,2047,3361,1952,5283,4046,1420,691,3999,2779,5216,1465,2702,3720,1904,2429,1419,2070,4944,3458,3491,5531,3249,1396,2439,5398,3285,5234,997,3353,5233,5124,396,2424,5589,2380,5706,1786,2451,3633,3632,5555,4884,421,1613,5126,707,3716,5411,4211,3241,3193,3494,4086,1736,5488,2412,1025,2150,1846,2128,2059,5603,5855,3192,5036,1224,1230,1228,2425,5281,5037,3258,3194,838,2788,1119,2226,3762,5143,2286,1949,3234,5539,3191,1620,245,5503,3177,4573,4051,1870,2445,3217,2116,2061,3244,3604,3490,414,4585,3801,5727,1897,5106,3603,3223,5389,992,5530,4034,2062,1844,1633,1490,1784,3942,1791,2317,3922,2153,4001,4002,5543,398,2797,2767,3189,2446,2761,5268,5278,1350,3554,1080,5113,244,406,2337,5558,3243,1645,3570,3077,3571,2753,2408,4050,3222,5617,708,3943,2301,4926,1610,3224,5887,3617,130,3648,2478,3647,1257,2223,1956,4676,2751,95,5368,4236,3125,5979,417,2130,2356,2737,2041,1189,2138,135,2379,433,1340,1190,1342,117,3894,5542,5902,5107,5210,5056,3257,3219,1253,87,1400,1399,5189,5404,256,122,3662,2171,1019,2298,3340,1890,231,2688,3212,1934,412,25,5511,5048,2085,4886,4271,5538,2231,2736,211,1623,5988,2399,2155,2454,3448,3543,359,3172,5064,465,1669,206,18,3291,2305,4006,1090,1462,121,3262,1662,484,1670,3725,1839,1485,709,5714,1557,3963,460,2716,3208,3230,5315,1245,2687,1699,2845,5100,2807,261,2487,3213,1611,2320,2763,3080,5599,3227,2717,3392,5940,6044,4040,2217,5917,2734,5534,4056,2176,3989,5860,4225,1787,119,216,413,118,5775,1094,2104,1896,1117,3289,2432,4681,4344,3268,3239,1095,2172,2386,3286,5957,2452,3835,3206,2300,4602,3347,1978,5052,116,3480,2282,811,144,1118,1589,5127,2453,638,5621,1384,2530,5774,5616,226,3660,1845,3225,5256,5815,2139,5854,4883,3246,1933,2723,5066,1162,145,1842,3840,4709,251,5566,5358,3921,2077,84,5407,3367,639,596,5406,3371,120,1950,2281,1831,1243,1908,1390,13,3124,2377,5586,4349,2586,2278,3229,5403,265,1918,3293,3303,5821,1796,2136,1389,5665,5861,3613,2163,690,1083,5387,1154,191,4224,5391,4736,5192,5169,5168,3542,4020,411,3881,2406,1625,3602,115,5133,5642,96,5386,5142,2809,4005,4651,4037,3587,837,710,2698,920,3715,5707,3331,386,246,5397,5554,2131,2319,5246,6018,3263,812,3362,3569,3799,5244,5367,1576,805,6043,3674,5708,368,6007,3041,3187,4025,2185,5427,3473,485,3079,2893,4014,1826,5054,5593,3180,5790,5084,3356,3904,1267,5190,3185,5777,1360,1385,5592,5596,2953,3297,1646,2180,2642,5959,4570,3426,4202,5170,2983,5264,3553,1086,3619,5101,5430,312,5431,2994,3580,4036,1642,4763,2357,3714,3309,1489,2411,3620,3771,2069,2120,3531,4054,3150,1258,2238,3149,6031,207,2308,5823,1156,5383,1391,3156,806,2123,5134,5108,5114,3038,2697,4881,3163,240,5365,4649,3143,3364,2467,5858,2117,3322,4024,5773,5408,1044,2258,3254,2419,1924,1928,681,2082,2348,5263,5549,594,461,239,3337,2169,404,5297,595,1155,3586,5212,4015,5886,4300,949,1783,5684,4817,1643,2972,1046,1049,3056,221,640,1040,2747,257,131,3728,3868,3867,1330,1577,405,545,5915,5918,2576,1626,5128,1898,391,388,4601,5082,4879,3181,1346,315,919,1636,5405,227,4057,3176,5418,232,252,97,4896,187,3174,5550,147,941,3988,2400,685,5122,1045,2206,3316,3278,3349,2539,1616,2026,172,3916,3745,3744,142,5551,146,143,3798,114,5591,5667,1024,3889,2090,1472,4870,16,2637,1656,5033,2388,5651,680,5464,32,6017,3532,5636,5567,2724,23,212,818,5788,104,4615,247,1240,1476,134,3306,1122,1401,466,3238,5755,167,1866,17,3016,1339,2414,3113,4351,4948,5243,4705,5409,4618,2486,186,168,3941,178,4569,1128,635,5241,679,5251,636,3944,5515,3678,2066,2181,4013,1621,155,262,169,300,1054,185,2201,1591,2368,5888,3269,4358,4873,2950,1121,6005,192,4007,362,2328,3424,170,1193,156,176,181,2651,301,5662,3266,2746,4258,5364,2612,3393,5424,5250,208,5585,513,5818,3949,486,544,3343,5239,514,302,2581,1858,217,1856,184,5729,3284,2314,2525,1127,836,2484,4568,3146,5863,3252,258,5141,4619,1923,1194,2617,182,157,2249,3078,175,2894,233,2776,89,2175,4905,4613,4708,173,4125,1020,4949,228,2347,5402,464,88,5698,3359,5380,2054,5421,4813,3906,462,512,615,1017,4935,3248,616,5721,467,3112,5396,3918,5075,3850,3781,1917,2520,248,2466,4063,113,303,3334,2595,1795,1047,4121,160,397,161,2750,5292,5757,222,4880,299,298,3993,4895,5722,103,3071,3090,1120,3927,1192,6032,1261,213,5361,1081,1475,593,3456,5781,5294,3863,3862,263,158,177,4029,3321,3896,3717,98,1615,234,1359,3940,2058,3830,4704,3279,3374,2757,183,904,5615,159,241,5073,617,689,5797,2632,2835,3956,253,3345,105,193,2846,5266,2318,5045,4680,2071,4906,5514,4654,5862,5852,4053,1601,4299,1268,3962,3957,1631,1254,5399,5377,4092,4294,3829,1641,5768,5581,5355,3891,218,3422,5118,2203,532,4058,677,3363,34,2810,1945,5653,5121,3438,3795,55,1677,5482,3793,4900,2354,4226,4023,419,2152,2961,2398,5635,3002,835,4144,166,5345,4117,1565,5680,3144,4732,3198,3642,3499,4735,4035,3201,4816,6020,4812,2205,3282,5237,1082,3310,163,4572,1757,4600,489,5633,416,592,5723,14,5138,3166,3416,5112,235,180,678,1830,4353,5117,266,468,4887,5764,2340,3712,637,3978,608,3169,1000,311,3384,296,164,3952,2074,3365,1755,3783,2168,1861,141,2556,5681,1566,5812,4580,4907,4309,3919,4624,676,2079,1477,148,5137,292,297,4626,4067,995,3888,316,591,2505,5065,42,531,1005,2561,2204,5182,242,2299,3640,5825,2494,4979,140,174,1323,3049,290,3638,112,223,1388,5111,102,1588,3261,3242,1549,2500,1488,294,5720,1609,1748,291,4963,2756,165,6027,3846,317,4611,5035,2190,194,1091,3330,4653,3489,1619,4731,171,1627,3971,4703,5443,4874,515,2896,106,3815,3287,2815,3634,5454,1843,5442,3879,1329,6033,5453,1824,1164,5254,418,133,4579,4762,179,1361,3814,4227,2653,5664,2485,4293,5747,1246,5004,2760,3828,5487,4113,5688,293,1895,5992,3631,607,295,1785,3834,4989,1629,684,151,488,1903,5822,1048,2360,3951,1351,4008,1191,2785,5963,2378,2803,5034,4988,5374,1751,243,490,3672,3173,1705,4707,4596,3255,1902,2710,4759,6034,5644,3419,2338,289,2184,154,1085,1703,3657,150,612,4076,36,5605,5938,152,4577,4586,4965,3641,491,2009,2279,5685,153,5832,3366,4337,606,1328,5490,3313,162,3975,5907,5602,3803,1655,4811,2627,1614,5639,2160,1600,2968,5608,5289,1158,2701,3920,5600,5610,4964,2179,3875,267,195,3972,149,5656,99,3298,111,5654,807,2341,633,4696,2257,5819,1527,2802,4285,4815,43,3752,3773,4790,809,810,4228,1012,1560,3304,3724,3897,2361,943,4929,1740,5943,2330,1050,5341,492,6035,915,2200,4758,2413,403,3339,1829,4947,4627,5641,3766,2468,1889,2221,5296,2037,2144,6038,2035,2670,1548,1394,2256,2571,3899,2978,3270,634,3418,1628,2390,590,4266,4591,5971,5927,1534,3065,236,4016,3630,3108,4297,318,3659,3983,288,4610,3739,3740,1648,6011,1827,1862,940,4298,2290,1123,2990,268,613,110,1828,4238,5606,1003,1030,3907,1841,4734,898,654,3022,3909,6036,3447,951,1815,1961,3658,38,5038,2691,4622,683,3394,2118,5013,4588,2794,3000,4307,2358,891,4878,5746,6028,3950,4978,510,4592,3296,4614,3816,3164,3800,3807,4730,589,4578,1953,2722,2405,3121,2044,3264,2541,6004,3066,4030,688,107,3521,5526,5228,6037,2515,5132,2187,4354,5524,1762,3076,4128,1307,1654,1201,196,4109,1706,5624,319,2844,5014,3994,6026,1759,62,5003,2469,2255,3072,632,5994,3845,6024,890,998,3806,605,3608,4126,2847,4950,1526,4235,2029,4761,4608,3402,2470,3500,2365,4933,942,3421,4247,4363,1964,2042,2353,3765,2890,3584,2049,3503,5062,101,5131,4648,5990,832,614,5891,5894,109,4064,1552,4026,4230,5051,3235,3067,2796,3290,3573,699,3768,993,2036,2681,2028,3844,3188,880,4688,74,3087,3086,108,1758,1754,4700,1756,4311,831,542,1412,1594,278,2040,1638,6002,4231,918,2385,4621,5660,5895,3183,3329,1530,892,269,277,923,237,4122,4804,5998,1644,5646,530,2038,5287,287,4320,2065,5965,1812,1624,5753,100,2886,6039,5017,279,3215,2084,73,5658,924,913,54,682,5609,3609,1640,1700,914,1929,543,3033,188,83,400,5193,3926,4248,2479,1602,1200,3336,1313,5898,2537,1930,2597,911,3738,2471,1814,4757,5866,906,1587,2310,4268,3520,5892,3692,1782,3054,5546,2177,4239,2210,3506,4059,1312,5839,3414,3796,4118,659,1744,3420,1607,4669,944,5679,197,2121,5116,2195,3824,1847,894,2212,320,2472,45,4364,4234,4232,4699,5911,3178,2097,2183,6009,1084,845,4810,3210,5750,889,4702,5893,2849,4372,1901,5282,270,5342,59,3005,929,847,2917,3707,597,664,4028,1743,3088,3063,4077,1022,651,3294,5904,1139,2199,658,780,5291,326,907,3081,3307,1886,4661,238,4369,1780,3273,1867,3556,3555,660,611,2672,4124,3074,511,3544,3545,63,5999,1124,844,833,652,3691,819,3964,846,3060,4687,3823,2805,1750,3004,1920,1039,2125,2209,4777,2196,2043,3684,1218,5693,4576,5683,3050,3052,1800,3046,905,5787,1411,610,4365,4341,2825,3432,4317,3976,3923,508,526,474,4259,3317,1018,2270,3861,2094,1813,4336,4642,5813,3711,912,675,3945,1408,6022,1778,2420,4597,5754,5506,2657,895,843,3408,64,3635,58,2259,4571,3537,2448,3026,3171,6041,896,897,3538,5725,4066,4114,4221,1946,2982,2960,1833,1428,4807,5985,271,5903,1544,1234,524,4675,3649,286,357,4951,4872,3048,3925,3644,3643,5187,2592,2473,2441,4567,2050,3053,2474,2480,893,822,4750,3866,327,5941,4808,1087,4634,60,2475,1129,4789,1092,2842,604,927,2863,1592,4368,899,4633,2045,3706,1781,5061,1459,3661,5997,5996,2888,1027,2434,3073,5986,2671,2854,1337,1887,4081,4888,1776,1522,2051,3083,198,507,4660,1914,1770,5784,473,355,2958,4805,1456,3606,1801,2601,553,276,3756,3040,280,5005,53,5669,4806,1864,5912,2679,4314,5946,2427,804,4672,2476,3024,285,2837,1763,328,3936,1761,321,525,1103,4729,3299,3682,4223,5063,1779,5678,3743,476,4956,5525,3084,808,3283,4936,4788,1652,1936,5338,5831,534,477,1657,4583,2052,2536,1959,1056,4686,2986,5763,2101,1729,2889,272,3042,1055,2945,401,358,480,4673,2216,2593,446,3639,529,1802,1455,1894,2909,2103,631,5652,1869,1348,5910,552,1749,2649,325,650,4595,3937,4584,2925,4999,4240,4315,3161,4105,399,4769,948,1029,1429,2423,3109,2208,3817,4780,4018,281,6015,2980,472,5024,554,5255,5991,1792,4110,3804,4068,3977,2923,3019,3075,1381,4980,3974,3973,478,4992,5744,3946,630,2939,3774,1478,3596,588,3588,3589,509,1977,1838,1940,3271,2410,3089,4646,2738,506,665,4786,402,4756,5697,4966,663,5044,5696,6040,4742,72,479,3784,2869,61,4251,5805,2853,284,3082,1014,1159,5830,4250,2262,2188,4753,283,4575,1427,3653,1639,1806,3651,925,322,4723,4038,5796,3726,3327,44,1433,3685,4625,82,1138,4069,339,2689,2924,471,834,1402,6008,3652,4120,1967,5185,3061,4184,2189,5284,4679,5492,1760,4273,324,1954,4659,1457,878,1597,2778,5905,1207,1596,830,4356,5901,1747,3062,1742,821,816,1493,4781,5074,4809,4785,5339,4901,5692,4685,1860,273,1347,5527,2868,3320,1697,2956,4658,4620,2725,4582,2364,2740,3776,851,5982,674,2402,338,1974,5021,5030,1430,2325,1922,4889,5968,5964,2126,282,3854,3855,1891,2404,2852,3280,275,1820,4796,1015,5091,4200,888,3598,3301,3758,2148,2648,2193,4004,609,5346,853,523,5032,3032,4631,1723,337,1832,3705,2954,1441,1051,2860,5614,4943,3965,2989,329,2124,2243,936,1799,1938,4220,3498,2834,5699,4303,4218,2113,503,828,1966,4928,3051,857,528,4902,3122,4241,1836,3786,2129,2133,1653,1424,3646,5791,5136,4701,887,1907,2708,4877,4996,3011,340,3104,603,65,5916,602,360,330,2266,5020,559,4017,323,5071,1888,5028,3057,4754,4726,2841,2952,4312,5119,4768,902,5029,5989,2996,1944,2967,3645,5686,341,2707,336,1491,4715,4009,4779,2813,4127,1733,3162,2367,5970,4027,4242,4767,4632,5015,4645,946,1741,3650,470,1498,1650,825,4741,2291,4264,3565,3012,928,1583,1837,827,52,2096,1066,2214,4740,2241,5648,3597,5638,3593,2851,3552,4814,1292,4123,522,882,304,6010,3683,826,3797,3592,1304,852,1948,2246,4201,2384,5139,3576,353,516,938,5188,1023,4995,2081,4898,2137,1454,2285,2127,4350,1410,2908,4728,3858,856,1727,4003,3110,4306,5981,3526,3527,2545,2820,1335,1554,3564,4706,4727,2699,2046,885,5931,4106,1426,274,5070,1065,6006,4990,1297,3123,1717,1038,1298,4355,2119,2912,1973,2265,823,4587,501,657,541,3064,1208,335,2146,2132,5770,5027,349,2715,1708,1409,3948,3549,4698,3969,4181,2814,4305,4674,4253,3010,2944,3043,334,3679,947,350,1028,4558,3548,2870,502,861,4652,3575,886,1900,550,4714,3676,343,4647,4119,2313,4875,5782,4308,49,5906,4733,5792,4794,3446,653,2075,3541,3938,342,1439,2393,686,2107,1593,2787,4367,5899,1425,66,3018,858,332,2864,2269,71,3159,815,5547,1483,1752,4188,1957,622,2373,3509,2988,1320,3693,666,1958,4795,3020,2294,4752,533,5995,2333,347,1899,345,859,2089,2786,75,5310,4891,85,2248,4116,4599,5031,704,348,4904,1947,5343,2748,570,351,5135,783,673,5092,2966,4428,2946,3560,1835,333,4739,331,4713,1555,346,3559,2105,5914,2274,5846,567,5924,753,1731,1582,850,814,81,5344,2198,1840,4784,4115,1137,500,5634,2370,549,48,1386,3563,4755,1637,5668,2955,68,2974,2240,4252,4783,662,1470,1275,3025,1440,5751,1458,539,1572,4617,4712,4204,587,535,3035,5919,1962,4694,5980,1160,4598,56,551,2981,5820,1746,3027,3686,2867,5762,2215,3111,855,2091,352,2192,3003,344,4071,46,1575,2087,4760,3664,4244,5583,4000,785,661,4010,4361,801,2268,76,4078,4766,5745,2261,1963,930,5913,4319,5089,6013,3722,1321,1195,4112,4691,1725,2220,1651,2680,5827,6042,5817,817,3663,4934,1093,5184,1037,69,2088,5504,571,2421,4850,2926,4108,1819,883,1863,536,1798,1453,3998,1581,1395,5340,5286,2887,4295,782,2840,2039,4823,1432,4799,4482,2709,854,1818,797,2098,3694,2073,3085,2833,1906,5952,1816,2134,2304,5974,4640,527,4644,4946,2963,520,3735,1883,2100,4563,2959,57,4342,4313,558,4775,2032,800,505,1634,2254,848,4667,70,5311,568,4849,4243,2053,521,4101,2197,2033,3117,1147,2250,4111,4820,781,4962,3539,5093,2156,5088,5958,4831,1274,3047,4594,4802,4725,4822,629,2911,3036,3098,1088,1336,5944,5785,1993,3044,1885,4096,1135,4210,2839,3068,3677,4876,2002,6003,4213,47,1262,4748,499,3021,2030,4509,5840,1349,2005,1972,2218,4590,4107,566,577,3522,3523,5956,4955,1564,1916,1585,5060,5305,2831,5637,77,2758,497,15,3721,4347,1406,4100,80,1210,4986,2832,4616,4152,2976,4182,4800,4148,4821,2998,4829,5072,1379,4401,1114,585,1217,931,1474,4190,4070,435,4499,4566,5505,576,2164,4984,3014,2971,498,1537,4991,4104,2252,1994,1032,922,2948,4357,2862,1618,569,4844,1405,1438,4637,849,4684,4256,2823,3594,2271,4156,3871,860,4192,5507,4255,2941,4778,2135,4589,4376,4362,4360,4185,3029,4371,4671,1442,4987,1586,586,2008,6001,6000,2178,4222,1932,1605,647,4835,4203,932,1034,4695,4197,649,67,1561,4692,4834,2350,2947,5752,1969,5018,778,5115,1052,4953,1368,4993,2970,725,572,5926,1989,1211,4160,5083,2191,1371,1216,573,2977,1608,4958,4867,4102,1089,1777,3805,5026,1915,4176,5140,457,2083,687,4180,1538,78,1104,5016,2933,2001,5923,2937,2934,2086,5006,2836,1480,1805,2344,5928,2881,1206,2916,5080,4721,4472,789,4455,1968,4164,3037,900,648,79,4172,1991,4168,4079,4925,5306,5012,3105,4793,4229,2914,2949,5900,1486,2856,2999,4194,4536,4803,4348,1715,2816,829,1533,2283,2342,4075,463,5327,3550,1108,2920,799,786,4630,2006,4206,671,3015,1721,2804,4945,4427,5816,1985,1804,5795,4343,5298,2910,1112,5011,3034,5007,619,2690,4052,3013,4191,5329,4073,575,1990,5966,1997,4195,2884,5050,2992,1435,5328,3811,1574,361,5325,2993,2826,4882,901,1334,2731,1531,2093,842,547,1294,1125,1064,779,3534,3533,1872,4638,2014,4981,4845,1719,4657,5009,1584,1892,1745,5834,1072,1753,1113,1649,3870,5042,3591,2806,5759,574,4848,5001,5749,1443,3007,1541,1998,4335,1387,1986,644,4553,1150,4209,1136,4301,3547,2777,4833,4340,5993,2985,2245,2263,5043,584,4183,4985,4751,3700,51,1487,1098,4776,487,1145,862,1074,1247,4689,4581,672,2969,5090,1931,2755,5079,2407,1071,2824,2818,1532,3091,4103,1077,3747,4847,4498,4982,5975,4960,4212,4772,1324,4664,583,5977,5039,796,2700,4923,4198,1704,3928,1222,1431,5896,2416,1325,1134,2735,1102,726,3808,519,720,518,4471,841,2919,1148,1149,4897,4508,4526,2706,4334,4523,5925,1076,6021,2542,3701,1144,2588,4683,4913,1868,1033,2213,2665,775,495,4749,3028,2838,3558,5069,494,777,582,4328,2067,1008,1482,4207,1061,1075,2795,2021,2861,4668,1110,4216,546,5293,2034,4636,1380,2219,4504,2609,787,1111,2022,4074,548,656,1221,4279,1109,1524,1035,2076,4641,3116,4843,1473,4418,4039,3126,1713,4856,4787,4932,4481,4961,2025,2975,1999,4496,1327,5312,5303,4142,1987,2745,776,4143,4791,5326,3982,3981,1711,4920,2885,955,4217,4497,4738,1434,5308,5647,5010,5978,4663,2194,4346,1146,1996,4764,517,4912,4711,4290,3748,4745,3561,695,790,4629,1220,945,1540,1277,4797,788,4765,2102,4552,1416,3734,565,2859,4193,2661,917,2931,2247,3813,4819,2811,4773,2991,2964,2918,496,643,5640,2382,2639,4196,4665,1305,4470,3045,4524,3006,5120,1995,1309,4550,1469,5969,2987,1468,5307,4535,4662,5826,2906,3423,2855,4782,2812,5984,4860,504,2857,1769,2020,2940,4937,2817,1913,5330,2650,4151,5743,1067,621,4469,1151,1418,469,1223,784,1803,4840,564,422,3812,2801,4818,5929,1235,1789,628,4746,1209,1982,1316,2883,1315,1546,5677,711,1983,4415,3120,1595,2324,4737,537,2605,4099,623,2843,5081,540,2000,5643,4656,3535,1317,5922,1308,2538,2004,4561,4186,6023,4179,2819,3540,2003,1301,5748,2930,4562,2965,748,933,3723,3536,3160,5908,1471,1525,4219,5811,4090,600,718,2362,3688,4858,1536,4175,4635,4417,4339,4400,4072,4770,1105,4690,1073,4416,1797,3023,1570,1196,1289,4792,4370,5909,734,2608,2092,5008,4548,2928,4155,3009,3732,3731,1547,4525,2749,1988,6019,1370,4908,1007,2598,1417,2099,1445,4551,3417,4278,4147,4167,1326,953,5594,4940,4534,4159,5258,655,4282,3873,3872,556,1821,4559,4171,4846,2927,752,4565,1460,4507,1451,2418,5300,4163,867,4555,2578,746,4842,4150,4098,4454,3546,4710,2251,4146,2858,4718,1383,4722,3785,3127,4199,560,4095,1481,3953,5776,5961,1115,2899,2711,4091,2644,4154,627,2871,5167,4719,4915,3101,626,4771,4743,4423,2903,4375,3775,4697,3595,2957,1563,4208,555,2553,4864,1817,4374,3930,5976,4094,719,903,4522,4478,4477,2759,4274,4909,884,4511,3590,1716,5930,4744,4538,3513,747,4830,1775,1100,2068,774,2874,1290,5313,5087,3690,2007,4205,921,1415,1449,2962,4270,3767,4158,4388,1285,1036,1573,3411,4521,4480,4510,4430,4333,2664,4899,1768,641,771,4453,5632,5960,3689,4474,2821,4284,2233,4463,1357,5921,578,2532,4564,2904,4930,2729,2583,2659,3832,4215,4938,4097,3551,4166,4852,4162,538,4426,2915,4137,2594,4983,4174,791,2730,2935,3857,4170,1376,4473,876,2186,5301,1452,4549,1505,1345,2936,4283,581,4918,3518,1423,3761,3008,4724,4484,714,1358,2739,1682,2875,4178,3709,1269,1807,2902,1450,1099,1078,1976,1101,768,3877,3878,2017,1314,2024,5758,4853,1437,1393,1011,4483,2095,743,4839,1319,4442,5023,4457,1484,2012,1811,2106,2016,2879,2522,4445,5973,4832,5295,1279,1422,2019,2552,2865,601,4546,1771,2921,5000,1272,1730,4967,1808,2872,2873,729,1444,793,1822,4391,10,3524,4798,2614,2634,3001,4277,2549,1366,1270,2880,3727,4593,1062,1421,3703,4399,1728,3517,1447,4494,4467,3429,4412,5954,4717,4939,2876,5078,1063,2629,5059,2013,3929,2878,3031,4994,2822,2984,875,4500,4931,5040,4403,4836,1773,4475,1606,4276,3528,4854,2640,2877,1543,2619,4390,1364,4827,1310,4505,2932,4976,2979,5019,3562,1774,4670,1553,4318,792,3529,1079,1276,2938,2558,4065,4537,5025,4857,4456,1448,5920,1521,2866,4280,2563,1302,979,1273,3530,580,2011,5628,1271,3557,1992,1446,1356,4554,1229,1718,1523,3831,2943,2244,1520,730,1503,754,744,4476,4429,1286,715,4444,618,4450,1701,4387,4424,4544,1496,3939,733,2527,4862,1280,4272,4458,4916,798,4459,1374,5807,4998,2544,1975,4825,1299,703,795,5739,5290,1884,4532,1060,1511,2891,2507,763,2654,5673,1495,1542,1970,4528,4861,697,698,696,1278,1684,625,1291,2023,4214,727,1212,5972,563,1318,700,6014,762,4826,1097,5288,1300,1539,1140,769,1698,1504,717,4716,712,1857,716,1955,4443,4921,2827,4407,4492,3935,2015,4381,1732,4389,1009,669,4408,1545,732,4541,4402,1693,1562,4292,1971,5302,2997,1367,877,1288,4495,1068,2584,1506,5949,702,1550,4527,2018,4493,5772,1303,934,4539,978,4973,1306,1764,4489,4431,2417,3512,3511,1282,5948,2010,713,4377,4501,4917,772,4420,1106,646,1355,1695,1287,1724,4506,4643,1377,1436,4404,4304,11,1767,4419,1519,4859,2922,2603,4468,1479,4275,1284,4291,1293,1281,3810,956,4485,1070,2559,3114,1168,41,4093,2564,1722,1265,1960,3933,4639,1283,4451,760,6012,4488,493,4487,2768,758,645,4135,4383,4486,6016,4560,4289,1772,794,937,4380,2502,4386,4392,4330,4161,4157,3934,2897,1965,4149,2898,562,2907,2905,4841,2630,4145,1502,3477,1529,2508,4345,4153,770,2900,4503,4169,4177,4173,4165,4801,4530,1603,5897,728,701,803,4080,5068,742,4529,4517,4269,2620,5285,3708,4421,2913,4479,1322,5049,1250,4914,3190,4502,5587,4393,670,1694,1187,2547,3510,761,722,4666,1692,4439,8,12,2830,4837,4440,4331,4452,2575,3,579,757,4140,4465,1042,4410,4460,5590,4540,1516,2829,4409,731,741,667,3525,4838,1010,4824,4435,4974,7,2503,5670,2615,4491,2611,1143,3980,4385,561,4557,668,2892,1765,4359,2585,4287,4464,1604,755,2901,4851,4288,4434,9,4187,1921,2392,773,4972,4462,1375,4447,1169,4397,2607,4405,557,3687,5231,4512,4519,4545,4461,4395,1225,4446,4466,4413,4547,4425,2658,4139,1373,1528,4411,4141,3516,2567,4556,4422,4396,4922,1691,751,4919,4134,3515,2555,977,3486,871,6,4774,1174,4366,4490,2791,4441,4910,2929,2566,1372,4189,4543,1295,3519,745,1709,4969,2606,4513,1726,1905,1766,2510,4398,4952,475,2882,4384,50,1178,3469,3462,1509,5737,4720,642,738,2548,4514,2517,2349,1535,1392,1369,4129,1559,5304,1937,1215,4136,4378,723,3802,4406,5,1518,2369,4138,4322,1219,2782,1712,750,5099,5577,1984,740,2645,2573,2942,4911,1311,4828,735,2533,749,2668,1514,4868,3754,3755,4533,1515,4,1,2828,1710,4286,2622,5579,1142,4449,2589,5002,1133,4379,2516,3986,3987,4323,3985,4382,2554,1707,4352,2543,4131,1173,4281,3876,5625,1205,2570,4520,756,2765,2528,2663,2572,1492,5022,2322,2579,4394,4448,3992,3991,4437,2582,1363,4977,3496,3954,2666,4438,1238,2676,2562,1494,2372,4866,4997,5098,2662,2560,5058,2551,4515,1714,1690,4089,2587,1107,2577,5262,2312,481,1497,2590,2568,935,1213,3809,1197,1183,5208,2293,767,4518,1069,1499,1720,4531,2253,3030,2646,2574,1365,5077,5041,2273,3997,3996,2686,2529,1180,4516,1172,2521,4326,2600,1182,2638,2569,4324,957,2565,2511,4338,2655,4133,2696,2599,5299,4414,2635,2523,4436,4941,4954,2604,1378,2667,1177,2550,2602,5076,4321,1508,1517,2557,739,4433,1500,4970,2623,1141,766,2506,4975,759,2531,5987,2624,624,4432,2773,764,2546,1551,4325,4542,1512,3995,5097,1296,2534,2643,2513,4327,3681,2628,1893,1181,2514,737,1501,1507,1510,4863,1227,4968,3514,2660,3990,3931,2518,4747,1513,724,3092,1184,2519,4924,2060,721,1226,3750,3749,2512,3753,5309,1186,2501,2,2526,2409,4296,2524,765,2618,1171,3826,2610,2616,5095,5235,2509,4130,1170,4971,1214,4332,952,3932,5096,4693,2621,5094,5067,3924,2641,2626,4132,2504,4865,4373,736,620,1176,2580,5086,5204,2613,5842,2631,958,2636,2633,2625,802,1179,5085,1185,954,2656,4855,2728,1175];
TOP_FEATURE_COUNT=200;
FSLIB_TOOLBOX_DIR='/Users/Zhao/Documents/MATLAB/Add-Ons/Toolboxes/Feature Selection Library/code/FSLib_v5.2_2017';


