from abstract_quantum_simulation import *

num_calls = []
max_calls = []
for H in range(10, 241, 10):
    print("\n\n#####H=", H)
    diff = math.sqrt(1/H)
    vals = []
    for center in range(2, 6, 1):
        means = [center/10 - diff/2, center/10 + diff/2]
        print("########\nevaluating for means:", means, "at H =", H, "center =", center)
        vals.append(BestArm(means, 0.1))
        print("vals for this:", vals)

    max_calls.append(max(vals))
    num_calls.append(vals)
    print("all results so far:", num_calls)
    print("max results so far:", max_calls)

# delta = 0.1:
# all results so far: [[601363326, 459237726, 465582915, 488918779], [1063259910, 1087919634, 959134704, 1063323534], [1201744374, 994933194, 970545684, 957283134], [1035736926, 1066609926, 1040635848, 1110847290], [899614830, 1555423506, 1703538228, 1064936580], [1628184960, 1526256840, 1413541752, 1503533217], [1631290872, 1602985026, 1455811062, 1646902338], [1718624796, 1516213590, 1434759288, 1445104059], [1654398402, 1534906824, 1687045428, 1561287243], [1521281832, 1678341264, 1687045428, 1701550452], [1627027620, 1759677504, 1609080447, 1698192678], [1622915622, 1766141190, 1360748667, 1701550452], [1622915622, 1547584784, 2545554090, 1411558070], [2625338142, 1547584784, 2545554090, 1414729301], [2860000860, 2550003210, 2399471472, 1411558070], [3152860530, 2976340368, 2852079350, 2666428373], [3016995122, 2823546296, 2762386616, 2893873181], [2890793072, 2823546296, 2456092856, 2744432789], [2875457612, 2812970612, 2683537664, 2741261558], [2849679146, 2820188522, 2534097272, 2744432789], [2890793072, 3123603776, 2834154752, 2544115145], [2741352680, 2666923256, 2834154752, 2696661449], [2741352680, 2520840638, 2809767554, 2696693905], [2748570590, 2517482864, 2809718870, 2696758817]]
# max results so far: [601363326, 1087919634, 1201744374, 1110847290, 1703538228, 1628184960, 1646902338, 1718624796, 1687045428, 1701550452, 1759677504, 1766141190, 2545554090, 2625338142, 2860000860, 3152860530, 3016995122, 2890793072, 2875457612, 2849679146, 3123603776, 2834154752, 2809767554, 2809718870]

# delta = 0.01:
# all results so far: [[729233382, 548189904, 583601367, 592436599], [1739219922, 1730393970, 1574720952, 1846643154], [2165552274, 1618810242, 1628709612, 1568301063], [1668754818, 1742033658, 1669130472, 1903518090], [1436422818, 2722521714, 3198391152, 1701556560], [2962320210, 2679103344, 2395151832, 2650501245], [2980350960, 2930951472, 2593358742, 2984081085], [3063012240, 2666819694, 2417607528, 2273653428], [3023399658, 2685222246, 3178754160, 2948514144], [2515116426, 3013612752, 3170016690, 3191660358], [2942683218, 3265428330, 2933024934, 3195764304], [2951976498, 3351758010, 2354862738, 3195764304], [2951976498, 2713111664, 4801968804, 2396446481], [2951976498, 2709007718, 4797864858, 2396446481], [5728698048, 3037398224, 4298931264, 2396446481], [5978735264, 5945919528, 5486395508, 4934786651], [5811938606, 5442886082, 5205579442, 5699236451], [5524629086, 5442886082, 4517945096, 5196198911], [5524629086, 5442886082, 5282394896, 5196198911], [5515891616, 5438782136, 4779357356, 5196198911], [5524629086, 6113911466, 5454486686, 6385224713], [5178723122, 6118015412, 5454486686, 5282216213], [5021591546, 4589005772, 5486198996, 5282216213], [5021591546, 4589005772, 5430148328, 5282216213]]
# max results so far: [729233382, 1846643154, 2165552274, 1903518090, 3198391152, 2962320210, 2984081085, 3063012240, 3178754160, 3191660358, 3265428330, 3351758010, 4801968804, 4797864858, 5728698048, 5978735264, 5811938606, 5524629086, 5524629086, 5515891616, 6385224713, 6118015412, 5486198996, 5430148328]

# delta = 0.0316227766:
# all results so far: [[665298354, 491677350, 515390865, 540677689], [1146628500, 1175210706, 1032981732, 1139965608], [1289735532, 1072925622, 1064061822, 1030804704], [1126277466, 1149710376, 1118760480, 1192164954], [1011427056, 2091197814, 2387599194, 1149300474], [2016170460, 1907325564, 1790063904, 1881662703], [2167773204, 2136069840, 2032762656, 2332606887], [2244075924, 1913031120, 1809723936, 1824127641], [2044754490, 1913031120, 2369534298, 2085950442], [1898130702, 2216103420, 2361556608, 2373583332], [2014974570, 2305463160, 2139220401, 2385291882], [2166577314, 2453335044, 1728355128, 2385291882], [2158599624, 2082573428, 4088374623, 1787793142], [2158599624, 2082573428, 4095649800, 1787760686], [4817573160, 2381914868, 3611243280, 1791353687], [5069959334, 5057537496, 4897381178, 4224031127], [4909340330, 4573135070, 4475047940, 4969854131], [4643939670, 4569404210, 3669016400, 4485398927], [4627030094, 4565673350, 4414790720, 4485398927], [4635921410, 4565673350, 3930384200, 4485398927], [4643899100, 5219528750, 4580508740, 4759983125], [4320111584, 5219528750, 4580508740, 4395084749], [4167470270, 3734214590, 4609200830, 4395084749], [4159492580, 3734214590, 4702674200, 4395084749]]
# max results so far: [665298354, 1175210706, 1289735532, 1192164954, 2387599194, 2016170460, 2332606887, 2244075924, 2369534298, 2373583332, 2385291882, 2453335044, 4088374623, 4095649800, 4817573160, 5069959334, 4969854131, 4643939670, 4627030094, 4635921410, 5219528750, 5219528750, 4609200830, 4702674200]


