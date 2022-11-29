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
    for (int i = 0; i < cap; i++)
    {
        foo.push_back(i);
    }

    // vector<size_t> foo = {4, 7, 25, 20, 24, 16, 10, 9, 17, 8, 11, 5, 22, 23, 0, 2, 13, 21, 19, 1, 18, 3, 6, 15, 14, 12};

    if (random_)
    {
        // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        // shuffle(foo.begin(), foo.end(), std::default_random_engine(seed));
        // for (size_t i : foo)
        // {
        //     fprintf(stderr, "%zu,", i);
        // }
        // fprintf(stderr, "\n");

        foo = {777,857,177,711,315,632,297,64,547,824,924,843,909,302,517,327,351,464,471,407,757,253,405,986,893,252,604,161,102,505,598,774,66,215,929,625,593,704,754,167,343,504,722,412,296,631,121,826,78,617,980,440,304,179,335,819,907,742,545,890,958,3,319,149,133,183,976,337,560,69,42,139,173,901,591,616,77,561,381,978,137,310,275,822,866,70,802,531,72,682,261,334,401,263,103,863,552,86,933,842,700,54,669,184,192,191,409,135,707,860,587,160,496,499,28,362,329,308,251,365,357,807,247,465,710,101,152,479,376,235,460,112,569,692,393,923,483,340,872,656,309,728,331,193,959,574,37,745,668,197,259,352,73,404,347,13,358,370,88,400,731,712,395,36,610,301,402,790,435,686,355,939,998,285,226,590,706,287,536,943,108,640,418,449,48,142,469,652,917,388,265,750,119,621,172,645,477,240,423,732,186,602,232,854,793,228,918,359,330,635,416,363,993,538,200,781,224,442,808,307,953,109,169,932,782,605,262,949,979,0,58,948,825,196,527,935,291,972,881,486,683,867,399,684,888,719,839,791,468,210,558,519,266,132,654,439,851,671,378,786,834,856,555,63,163,708,280,419,250,679,755,911,462,688,493,487,31,841,94,403,427,450,387,900,481,816,780,190,290,458,743,211,649,467,114,20,125,264,433,577,181,914,578,294,787,549,736,208,445,59,973,209,944,626,217,117,22,32,273,767,29,830,630,156,835,71,779,157,534,877,945,55,612,989,556,551,920,971,513,24,83,38,138,926,349,984,906,815,609,954,212,606,886,847,887,49,255,868,611,245,375,608,735,434,85,798,348,80,582,516,751,738,568,205,239,812,463,410,586,730,500,236,624,776,23,520,188,799,637,21,565,642,507,87,553,389,341,201,783,146,2,62,699,473,981,800,763,968,832,713,511,414,965,512,934,122,116,543,258,332,579,805,134,941,246,44,120,653,508,233,869,279,603,339,836,544,148,896,904,698,619,397,19,541,764,472,485,634,758,451,206,443,234,314,67,428,922,650,489,813,424,542,34,272,749,716,563,369,219,131,510,916,628,845,213,810,417,51,659,846,39,655,806,202,41,964,974,377,940,638,454,127,350,837,267,599,583,864,829,550,809,658,539,249,955,436,159,674,715,691,76,89,748,827,997,687,260,151,383,60,106,380,657,878,322,651,336,903,105,618,372,475,573,346,446,185,667,277,797,514,622,643,600,661,594,81,129,633,784,164,147,766,214,371,453,79,136,110,98,985,905,873,950,874,882,93,821,853,268,174,796,803,515,522,421,680,811,636,488,804,47,231,501,557,494,52,957,947,831,155,56,726,166,257,589,195,130,739,756,396,852,858,982,221,554,503,175,128,689,576,491,970,384,170,820,963,379,4,289,238,353,912,977,913,662,408,765,46,537,15,540,425,670,588,495,572,898,861,850,737,394,615,672,306,928,666,529,641,833,761,999,8,566,204,326,770,398,75,902,723,823,849,894,885,530,695,524,162,899,441,466,168,293,53,97,995,768,773,762,994,769,300,476,391,752,229,222,660,884,237,862,33,470,298,361,760,373,794,45,227,220,859,915,744,18,575,891,567,292,30,100,367,991,11,333,158,368,714,90,946,639,506,844,919,930,456,171,646,562,241,360,218,581,321,942,35,875,299,614,406,952,386,338,701,801,382,564,685,295,759,281,57,354,457,927,270,696,366,996,595,969,392,459,311,775,282,437,17,342,892,725,248,9,676,426,243,271,702,910,118,818,286,727,814,502,126,548,312,284,724,855,570,571,675,43,283,532,199,778,225,68,936,324,694,962,104,623,303,344,690,718,356,546,753,16,145,461,153,140,6,74,444,7,509,429,154,644,498,521,14,931,987,871,65,26,772,25,956,627,992,620,95,345,165,198,180,951,525,897,817,40,318,320,364,492,420,448,729,187,115,733,455,746,478,203,390,484,717,665,960,99,720,497,703,411,84,983,216,879,415,741,143,244,789,937,474,328,230,92,681,278,207,747,840,966,990,176,883,607,431,374,788,663,721,697,876,317,82,242,678,50,848,1,447,592,967,178,785,648,559,323,961,673,452,693,256,438,705,523,925,988,276,223,771,123,597,288,480,124,111,150,613,325,795,865,938,921,734,144,880,535,482,141,61,889,107,189,194,792,5,528,96,274,647,305,975,526,91,113,182,828,596,490,709,10,432,430,533,313,601,385,585,580,269,895,12,870,316,740,629,27,664,413,908,254,518,677,838,422,584};

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

    for (size_t i = 0; i < cap; i++)
    {
        // fprintf(stderr, "inserting %zu\n", foo[i]);
        size_t clus = tree.insert(seqs[foo[i]], insertionList, foo[i]);
        // clusters[foo[i]] = tree.findAncestor(clus);
        clusters[foo[i]] = clus;
    }
    for (size_t i = 0; i < cap; i++)
    {
        clusters[foo[i]] = tree.parentLinks[clusters[foo[i]]];
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
    vector<pair<float, float>> centroids = generateCentroids(clusCount, 15);

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

    // compressClusterList(clusters);

    FILE *pFile;
    if (random_)
    {
        pFile = fopen("output-r.txt", "w");
    }
    else
    {
        pFile = fopen("output.txt", "w");
    }
    outputClusters(pFile, clusters);

    return 0;
}
