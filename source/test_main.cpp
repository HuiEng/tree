#include <random>
#include "test_main_cmdline.hpp"
#include "test_tree.hpp"
typedef test_tree tree_type;
typedef vector<tuple<size_t, tuple<size_t, size_t>>> output_type;
bool random_ = false;
size_t cap = 0;
size_t iteration = 0;
using namespace std;

static test_main_cmdline args; // Command line switches and arguments

void outputClusters(FILE *pFile, const vector<size_t> &clusters)
{
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
    }
}

void outputClusters(FILE *pFile, const output_type &clusters)
{
    fprintf(pFile, "seqID,cluster,ancestor,level\n");
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu,%llu,%llu\n",
                static_cast<unsigned long long>(sig),
                static_cast<unsigned long long>(get<0>(clusters[sig])),
                static_cast<unsigned long long>(get<0>(get<1>(clusters[sig]))),
                static_cast<unsigned long long>(get<1>(get<1>(clusters[sig]))));
    }
}

void compressClusterList(vector<size_t> &clusters)
{
    unordered_map<size_t, size_t> remap;
    for (size_t &clus : clusters)
    {
        if (remap.count(clus))
        {
            clus = remap[clus];
        }
        else
        {
            size_t newClus = remap.size();
            remap[clus] = newClus;
            clus = newClus;
        }
    }
    fprintf(stderr, "Output %zu clusters\n", remap.size());
}

vector<size_t> clusterSignatures(const vector<data_type> &seqs)
{
    size_t seqCount = seqs.size();
    // seqCount = 300;
    vector<size_t> clusters(seqCount);
    tree_type tree(partree_capacity);

    size_t firstNodes = 1;
    if (firstNodes > seqCount)
        firstNodes = seqCount;

    vector<size_t> insertionList; // potential nodes idx except root; root is always 0

    default_random_engine rng;

    // node 0 reserved for root, node 1 reserved for leaves idx
    for (size_t i = firstNodes; i < partree_capacity; i++)
    {
        insertionList.push_back(partree_capacity - i);
    }

    vector<size_t> foo;
    for (int i = 0; i < seqCount; i++)
    {
        foo.push_back(i);
    }

    // vector<size_t> foo = {4, 7, 25, 20, 24, 16, 10, 9, 17, 8, 11, 5, 22, 23, 0, 2, 13, 21, 19, 1, 18, 3, 6, 15, 14, 12};

    if (random_)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        shuffle(foo.begin(), foo.end(), std::default_random_engine(seed));
        for(size_t i:foo){
            fprintf(stderr,"%zu,", i);
        }
        fprintf(stderr,"\n");

        // foo={40,699,500,329,385,282,301,686,173,19,689,92,68,62,880,354,845,830,10,802,843,303,411,88,285,987,984,765,316,221,911,488,901,951,435,181,980,112,356,332,894,771,421,692,615,463,318,456,497,288,371,638,741,602,151,271,471,869,548,432,108,390,454,217,925,199,299,256,193,807,685,375,324,844,487,594,235,695,417,478,330,484,272,543,198,975,204,472,998,947,110,270,605,410,721,26,11,940,828,229,889,743,956,441,882,132,514,725,854,450,722,683,20,923,508,915,351,922,876,237,575,475,759,39,333,24,84,31,789,823,750,215,423,655,718,73,670,736,746,531,129,321,852,982,986,155,402,392,684,629,631,148,997,216,883,654,748,156,109,874,562,640,787,509,369,568,597,726,8,979,553,116,829,875,319,115,994,888,195,814,890,973,367,836,873,416,564,690,542,769,833,590,461,136,728,103,774,72,399,460,311,404,523,710,276,228,944,603,231,446,719,289,466,600,751,358,782,7,539,704,246,955,989,373,1,969,777,752,797,547,687,153,549,763,53,190,281,30,709,990,305,284,698,864,536,760,189,325,749,328,81,374,322,384,757,614,209,826,102,63,570,900,494,776,739,76,556,926,174,652,154,669,649,36,386,457,904,510,809,978,226,858,290,617,163,688,850,628,941,144,239,315,242,928,770,100,254,846,611,827,983,114,280,896,188,458,433,988,912,85,277,335,22,647,555,146,742,813,415,607,772,754,23,768,907,831,278,13,661,659,778,898,767,906,124,668,35,521,735,334,834,218,409,623,150,816,453,961,16,142,21,579,52,257,264,234,512,761,6,300,230,717,855,596,364,551,370,55,527,805,591,445,279,470,881,145,137,78,643,236,609,504,886,913,810,407,51,9,295,464,283,608,706,372,918,486,967,587,820,592,675,755,342,895,310,641,125,222,682,101,161,738,519,700,927,642,995,249,403,302,296,818,503,847,327,291,507,693,197,589,113,598,723,593,583,361,897,406,626,61,159,892,630,431,293,552,724,796,185,495,808,345,4,637,381,903,496,824,812,395,46,365,2,524,379,138,268,498,34,343,424,931,870,206,786,412,227,942,492,48,938,255,200,800,427,248,848,344,440,662,47,162,341,65,584,656,339,191,560,121,744,618,910,613,168,664,362,902,660,482,819,540,572,680,964,451,117,715,694,576,393,413,438,443,502,247,135,366,714,939,312,943,212,921,595,788,287,304,520,82,822,355,511,673,948,574,389,914,963,878,506,999,681,442,166,210,59,817,713,838,147,815,54,194,932,338,558,917,872,139,841,784,382,207,891,533,962,186,930,867,152,490,974,753,232,126,545,491,474,530,839,863,309,211,954,935,515,758,244,793,164,538,360,477,991,733,529,672,985,952,266,213,336,588,893,945,340,550,267,107,3,929,376,368,971,785,953,25,83,469,862,165,879,792,977,966,465,219,780,380,924,849,363,959,60,120,674,261,658,273,747,414,899,425,996,981,972,69,58,179,263,884,238,396,447,756,349,44,118,56,606,485,192,220,716,566,479,720,405,50,401,783,420,806,143,29,842,877,350,243,933,346,419,885,75,79,439,604,140,201,825,182,400,205,578,851,934,130,697,887,275,522,134,171,525,260,89,639,516,773,729,667,274,970,835,250,920,627,561,861,245,708,537,745,176,821,97,707,679,781,187,557,727,518,837,480,127,645,388,64,665,292,678,620,636,391,223,546,650,262,111,581,599,481,535,730,202,513,183,489,38,801,418,225,149,429,383,534,937,779,646,871,331,936,857,378,444,653,580,795,233,49,99,86,610,196,571,41,476,663,33,666,430,449,868,95,798,286,671,467,775,957,27,865,93,178,625,505,269,644,398,43,544,28,74,992,294,252,705,455,960,387,633,359,141,976,326,251,123,804,422,119,517,308,632,306,766,80,175,677,66,790,128,909,624,691,437,203,462,541,184,157,337,314,98,436,968,15,569,501,77,468,67,434,612,601,712,313,180,573,616,791,859,158,585,554,473,18,740,320,258,811,676,832,105,224,408,377,459,711,297,965,170,5,214,732,452,794,958,916,586,619,298,622,734,702,91,71,567,104,253,172,265,621,762,526,853,167,347,701,106,703,483,45,307,57,950,908,860,993,528,905,177,352,493,582,17,42,840,799,357,122,696,731,764,169,532,87,37,577,737,428,651,919,133,240,634,866,131,648,353,949,803,70,94,559,32,317,208,563,946,657,426,0,259,635,160,394,12,323,90,565,241,96,14,499,348,856,448,397};
        // vector<size_t> temp{212, 677, 185, 710, 785, 1035, 512, 29, 690, 594, 929, 1016, 801, 683, 1065, 211, 390, 789, 945, 408, 806, 687, 302, 736, 66, 260, 649, 507, 771, 100, 72, 417, 599, 1010, 981, 439, 563, 484, 291, 969, 11, 144, 328, 362, 530, 588, 906, 403, 540, 146, 832, 96, 351, 149, 319, 254, 104, 986, 971, 777, 944, 158, 639, 414, 421, 189, 889, 611, 52, 863, 276, 161, 850, 941, 853, 22, 962, 565, 486, 531, 591, 124, 127, 210, 247, 1036, 815, 516, 749, 514, 1089, 462, 696, 31, 511, 536, 650, 763, 398, 379, 520, 1019, 356, 823, 373, 1063, 384, 267, 557, 350, 263, 383, 664, 891, 141, 358, 446, 692, 852, 866, 36, 98, 34, 524, 570, 139, 103, 528, 663, 503, 290, 933, 521, 1072, 975, 1077, 74, 829, 908, 265, 502, 497, 1015, 54, 874, 917, 85, 821, 416, 655, 733, 376, 91, 345, 722, 883, 122, 916, 743, 938, 270, 946, 605, 830, 313, 134, 130, 812, 1062, 756, 553, 1004, 231, 355, 1055, 729, 459, 487, 781, 311, 208, 926, 495, 921, 934, 930, 1030, 651, 1040, 839, 457, 537, 6, 509, 741, 389, 949, 136, 714, 387, 1027, 489, 798, 175, 219, 242, 256, 545, 898, 413, 152, 207, 773, 300, 449, 407, 738, 108, 296, 593, 264, 21, 622, 162, 62, 887, 814, 105, 79, 787, 166, 816, 470, 1061, 835, 499, 236, 970, 222, 221, 595, 151, 586, 783, 884, 577, 834, 923, 20, 28, 619, 43, 548, 755, 1042, 450, 478, 485, 227, 55, 75, 615, 238, 112, 993, 107, 713, 431, 653, 673, 168, 602, 1082, 286, 732, 872, 145, 435, 253, 643, 1031, 361, 1068, 1014, 468, 542, 626, 498, 354, 1029, 634, 332, 697, 871, 120, 997, 967, 335, 724, 681, 609, 966, 1079, 440, 183, 106, 554, 578, 434, 566, 510, 747, 337, 341, 752, 474, 123, 942, 841, 186, 372, 70, 876, 125, 766, 39, 129, 192, 2, 775, 802, 523, 323, 196, 334, 1025, 411, 438, 625, 49, 339, 700, 527, 885, 19, 481, 188, 796, 1086, 117, 652, 38, 840, 367, 1000, 948, 243, 1081, 466, 678, 228, 460, 657, 564, 301, 271, 342, 92, 712, 76, 456, 644, 241, 904, 86, 237, 1032, 371, 689, 820, 453, 613, 282, 568, 277, 953, 892, 534, 1003, 479, 894, 709, 142, 828, 465, 734, 215, 95, 306, 826, 443, 27, 206, 680, 1083, 232, 508, 73, 786, 1091, 119, 138, 719, 285, 1070, 377, 620, 10, 193, 7, 546, 1018, 702, 667, 58, 541, 701, 89, 327, 817, 629, 261, 784, 116, 698, 896, 688, 860, 101, 396, 668, 441, 574, 225, 859, 430, 169, 632, 601, 424, 1041, 1007, 288, 148, 822, 472, 454, 137, 980, 529, 659, 195, 268, 297, 899, 312, 581, 4, 903, 60, 999, 1075, 458, 121, 153, 233, 589, 473, 746, 910, 165, 316, 513, 394, 381, 776, 645, 173, 1071, 562, 1005, 375, 71, 932, 555, 308, 150, 772, 155, 87, 401, 630, 331, 901, 349, 99, 199, 24, 272, 902, 410, 273, 329, 1006, 984, 654, 17, 220, 287, 600, 30, 596, 333, 269, 352, 661, 726, 936, 496, 190, 637, 476, 842, 167, 1090, 675, 404, 32, 1046, 5, 847, 492, 110, 867, 406, 543, 811, 631, 1023, 201, 244, 303, 3, 275, 1059, 1066, 585, 788, 880, 77, 990, 234, 205, 760, 735, 204, 235, 131, 1020, 405, 873, 539, 592, 448, 567, 556, 289, 368, 1060, 51, 480, 587, 606, 935, 579, 769, 656, 610, 983, 365, 742, 177, 658, 865, 704, 624, 69, 895, 452, 477, 744, 59, 877, 132, 1056, 957, 279, 616, 44, 344, 793, 517, 538, 836, 685, 309, 745, 778, 346, 1052, 897, 451, 395, 1050, 703, 1044, 988, 1034, 156, 638, 409, 40, 179, 33, 283, 751, 397, 1087, 53, 184, 708, 194, 780, 42, 315, 768, 900, 45, 322, 647, 571, 419, 386, 790, 1053, 392, 956, 947, 359, 608, 711, 65, 258, 172, 176, 181, 582, 266, 964, 753, 483, 612, 378, 532, 280, 257, 500, 83, 813, 764, 730, 1008, 918, 230, 870, 809, 559, 939, 494, 427, 444, 1064, 490, 757, 914, 535, 0, 202, 774, 318, 958, 799, 111, 731, 491, 133, 135, 976, 363, 607, 353, 203, 695, 209, 180, 164, 1073, 1080, 366, 475, 447, 295, 940, 875, 795, 41, 1048, 982, 991, 170, 597, 576, 636, 862, 1021, 779, 428, 580, 1, 37, 791, 23, 825, 429, 627, 728, 996, 218, 293, 623, 81, 420, 864, 274, 515, 94, 1076, 759, 721, 461, 869, 82, 330, 455, 157, 590, 909, 522, 126, 223, 837, 748, 699, 937, 550, 357, 573, 818, 118, 423, 26, 886, 84, 12, 109, 720, 762, 868, 18, 213, 754, 171, 340, 525, 679, 326, 248, 469, 348, 977, 959, 1013, 694, 603, 229, 1085, 992, 1017, 952, 1038, 851, 246, 418, 501, 214, 46, 400, 114, 432, 666, 464, 385, 973, 640, 262, 672, 575, 633, 154, 831, 1037, 890, 968, 642, 504, 370, 216, 97, 393, 974, 665, 78, 187, 737, 670, 245, 493, 1047, 648, 1002, 9, 569, 907, 669, 995, 422, 950, 794, 255, 979, 128, 888, 433, 425, 294, 922, 198, 739, 314, 321, 347, 217, 572, 805, 765, 88, 800, 912, 90, 682, 163, 920, 707, 716, 833, 547, 845, 723, 61, 1069, 505, 924, 558, 115, 965, 93, 1074, 718, 925, 767, 843, 374, 1088, 64, 854, 304, 761, 552, 382, 770, 249, 808, 369, 25, 943, 662, 671, 415, 855, 913, 343, 758, 1001, 824, 67, 1009, 445, 684, 693, 931, 1024, 676, 8, 725, 893, 985, 13, 963, 191, 810, 846, 583, 519, 1051, 467, 686, 621, 844, 560, 226, 292, 338, 1045, 250, 298, 336, 614, 113, 827, 544, 402, 803, 1043, 972, 928, 882, 1057, 807, 200, 15, 240, 961, 388, 1026, 1084, 750, 955, 436, 50, 182, 437, 838, 1012, 1058, 617, 717, 533, 259, 251, 706, 646, 1067, 463, 878, 674, 989, 364, 252, 715, 63, 978, 960, 281, 391, 551, 56, 160, 147, 797, 14, 911, 488, 518, 857, 740, 380, 915, 305, 1078, 905, 80, 16, 35, 927, 482, 1011, 628, 224, 804, 299, 954, 727, 399, 604, 549, 998, 819, 284, 881, 660, 1039, 239, 951, 506, 278, 140, 197, 691, 174, 68, 47, 360, 858, 994, 1022, 471, 1049, 442, 48, 598, 792, 584, 412, 310, 861, 526, 102, 324, 143, 856, 919, 307, 705, 426, 178, 879, 1054, 320, 1028, 635, 782, 618, 1033, 561, 57, 317, 159, 849, 325, 848, 987};
        // foo = temp;
    }

    // for (size_t i = 0; i < cap; i++)
    // {
    //     for (size_t j = 0; j < cap; j++)
    //     {
    //         fprintf(stderr,"%zu,%zu,%.2f\n",i,j,calcDistance(seqs[foo[i]], seqs[foo[j]]));
    //     }
    // }

    clusters[foo[0]] = tree.first_insert(seqs[foo[0]], insertionList, foo[0]);
    for (size_t i = 1; i < cap; i++)
    {
        // fprintf(stderr, "inserting %zu\n", foo[i]);
        size_t clus = tree.insert(seqs[foo[i]], insertionList, foo[i]);
        // clusters[foo[i]] = tree.findAncestor(clus);
        clusters[foo[i]] = clus;
    }
    fprintf(stderr, "\n\n\n\n");

    for (size_t run = 0; run < iteration; run++)
    {
        tree.printTreeJson(stdout);
        fprintf(stderr, "Iteration %zu\n", run);
        tree.prepReinsert(tree.root);
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = tree.reinsert(seqs[foo[i]], foo[i]);

            fprintf(stderr, "\n found %zu at %zu\n", foo[i], clus);
            // clusters[foo[i]] = tree.findAncestor(clus);
            clusters[foo[i]] = clus;
        }
    }

    // tree.mergeChildren(tree.root, insertionList);

    // if (iteration == 0)
    // {
    //     tree.printTreeJson(stdout);
    //     tree.clearSeqId(insertionList[0]);

    //     for (size_t i = 0; i < cap; i++)
    //     {
    //         size_t clus = tree.search(seqs[foo[i]], foo[i]);
    //         tree.seqIDs[clus].push_back(foo[i]);

    //         fprintf(stderr, "\n found %zu at %zu\n", foo[i], clus);
    //         // clusters[foo[i]] = tree.findAncestor(clus);
    //         clusters[foo[i]] = clus;
    //     }
    // }

    // Recursively destroy all locks
    tree.destroyLocks();

    tree.printTreeJson(stdout);

    return clusters;
}

vector<pair<float, float>> generateCentroids(size_t count, size_t factor = 1, float mean = 0, float std = 1)
{
    vector<pair<float, float>> centroids;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    std::normal_distribution<double> distribution(0.0, std);

    for (int i = 0; i < count; ++i)
    {
        centroids.push_back(make_pair(distribution(generator) * factor, distribution(generator) * factor));
    }

    // for (int i = 0; i < count; ++i)
    // {
    //     std::cout << centroids[i].first <<","<< centroids[i].second << std::endl;
    // }
    return centroids;
}

vector<float> rNorm(size_t count, float mean = 0, float std = 1)
{
    vector<float> output;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    std::normal_distribution<double> distribution(mean, std);

    for (int i = 0; i < count; ++i)
    {
        output.push_back(distribution(generator));
    }

    return output;
}

vector<data_type> writeData(FILE *pFile, size_t clusCount, size_t clusSize)
{
    vector<pair<float, float>> data;

    fprintf(pFile, "%zu\n", clusCount);
    vector<pair<float, float>> centroids = generateCentroids(clusCount, 10);

    vector<float> x;
    vector<float> y;

    for (int i = 0; i < clusCount; ++i)
    {
        fprintf(pFile, "%f,%f\n", centroids[i].first, centroids[i].second);
    }

    for (int i = 0; i < clusCount; ++i)
    {
        x = rNorm(clusSize, centroids[i].first);
        y = rNorm(clusSize, centroids[i].second);

        for (int i = 0; i < clusSize; ++i)
        {
            fprintf(pFile, "%f,%f\n", x[i], y[i]);
            data.push_back(make_pair(x[i], y[i]));
        }
    }

    return data;
}

// read data, the first line is the number of clusters, the next few lines are the cluster centroids
vector<data_type> readData(const char *path)
{
    vector<data_type> data;
    FILE *fp = fopen(path, "r");
    if (!fp)
    {
        fprintf(stderr, "Failed to load %s\n", path);
        exit(1);
    }

    fprintf(stderr, "reading %s\n", path);
    size_t clusCount;
    fscanf(fp, "%zu", &clusCount); // check how many clusters

    for (size_t i = 0; i < clusCount; i++)
    {
        char fileBuf[100];
        fscanf(fp, "%s", fileBuf); // skip centroids
    }

    size_t i = 0;
    for (;;)
    {
        float x;
        float y;
        if (fscanf(fp, "%f,%f\n", &x, &y) < 1)
            break;

        data.push_back(make_pair(x, y));
        i++;
    }
    return data;
}

int test_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    if (args.random_arg)
    {
        random_ = true;
    }

    if (args.split_threshold_given)
    {
        split_threshold = args.split_threshold_arg;
    }

    if (args.stay_threshold_given)
    {
        stay_threshold = args.stay_threshold_arg;
    }

    fprintf(stderr, "split threshold: %.2f\n", split_threshold);
    fprintf(stderr, "stay threshold: %.2f\n", stay_threshold);

    if (args.minimiser_match_given)
    {
        minimiser_match_threshold = args.minimiser_match_arg;
        fprintf(stderr, "minimiser_match threshold: %zu\n", minimiser_match_threshold);
    }

    iteration = args.iteration_arg;

    // vector<double> seqs = readPartitionBF(inputFile);
    // fprintf(stderr, "Loaded %zu signatures...\n", seqs.size());

    // signatureSize = seqs[0][0].size();
    // fprintf(stderr, "Building Signature...\n");
    // default_random_engine rng;
    // // output_type clusters = clusterSignatures(seqs);
    // vector<size_t> clusters = clusterSignatures(seqs);

    // fprintf(stderr, "writing output...\n");

    // size_t firstindex = inputFile.find_last_of("/") + 1;
    // size_t lastindex = inputFile.find_last_of(".");
    // string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    // auto fileName = rawname + "-s" + to_string((int)(stay_threshold * 100)) + "-l" + to_string((int)(split_threshold * 100));
    // if (args.tag_given)
    // {
    //     fileName = fileName + "-" + args.tag_arg;
    // }
    // FILE *pFile = fopen((fileName + ".txt").c_str(), "w");
    // outputClusters(pFile, clusters);
    vector<data_type> seqs;
    if (args.input_given)
    {
        seqs = readData(args.input_arg);
        fprintf(stderr, "Loaded %zu signatures...\n", seqs.size());
    }
    else
    {
        FILE *pFile = fopen("data.txt", "w");
        seqs = writeData(pFile, 10, 100);
        fprintf(stderr, "Created %zu signatures...\n", seqs.size());
    }

    if (args.capacity_given)
    {
        cap = args.capacity_arg;
    }
    else
    {
        cap = seqs.size();
    }

    vector<size_t> clusters = clusterSignatures(seqs);
    FILE *pFile = fopen("output.txt", "w");
    outputClusters(pFile, clusters);

    return 0;
}
