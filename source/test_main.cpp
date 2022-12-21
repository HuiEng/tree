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

        foo = {526, 101, 386, 65, 741, 683, 632, 785, 765, 456, 493, 840, 181, 796, 956, 126, 497, 409, 368, 786, 4, 265, 234, 142, 612, 832, 99, 833, 394, 16, 871, 22, 278, 769, 999, 219, 728, 256, 731, 369, 344, 172, 2, 61, 834, 91, 251, 638, 984, 531, 923, 627, 527, 466, 163, 478, 859, 147, 32, 547, 652, 946, 319, 259, 366, 592, 541, 783, 735, 55, 922, 775, 97, 174, 414, 913, 803, 467, 890, 169, 364, 534, 248, 800, 795, 812, 161, 66, 374, 928, 845, 69, 400, 804, 407, 572, 606, 525, 872, 426, 938, 720, 390, 83, 819, 59, 582, 171, 518, 322, 336, 886, 944, 468, 491, 170, 510, 760, 477, 577, 18, 195, 301, 502, 546, 499, 481, 498, 296, 996, 935, 173, 54, 693, 979, 787, 694, 781, 780, 15, 666, 925, 184, 48, 470, 571, 215, 113, 711, 21, 435, 433, 628, 868, 197, 548, 564, 354, 911, 965, 658, 915, 324, 676, 144, 853, 216, 191, 901, 692, 563, 51, 62, 633, 927, 68, 975, 317, 106, 823, 160, 855, 185, 130, 236, 220, 727, 168, 749, 80, 427, 385, 970, 930, 496, 247, 784, 679, 186, 308, 544, 809, 92, 28, 841, 128, 441, 810, 203, 707, 507, 994, 47, 411, 193, 649, 948, 584, 703, 229, 939, 384, 408, 900, 560, 418, 119, 87, 931, 776, 188, 159, 918, 536, 487, 469, 177, 697, 532, 739, 597, 932, 280, 519, 774, 667, 131, 12, 653, 449, 545, 865, 877, 663, 183, 25, 461, 460, 851, 998, 713, 421, 677, 665, 298, 920, 669, 878, 520, 788, 79, 867, 72, 846, 723, 782, 231, 542, 554, 540, 822, 104, 882, 729, 207, 905, 480, 993, 399, 311, 348, 704, 211, 447, 306, 579, 725, 625, 835, 419, 568, 446, 926, 380, 642, 726, 450, 451, 327, 10, 600, 685, 60, 820, 114, 454, 74, 57, 108, 338, 23, 578, 289, 341, 405, 762, 750, 580, 132, 566, 807, 561, 488, 559, 141, 959, 537, 439, 969, 316, 5, 576, 759, 389, 958, 381, 373, 992, 618, 940, 549, 133, 989, 982, 594, 960, 631, 436, 885, 864, 93, 672, 724, 856, 543, 670, 38, 621, 770, 282, 453, 893, 382, 109, 955, 718, 402, 570, 242, 179, 300, 805, 180, 586, 506, 876, 392, 986, 90, 312, 270, 281, 388, 674, 883, 274, 978, 416, 777, 71, 910, 277, 966, 888, 224, 639, 689, 737, 997, 828, 722, 35, 852, 556, 607, 339, 258, 148, 315, 417, 415, 33, 830, 430, 401, 201, 39, 476, 291, 455, 736, 264, 943, 591, 398, 58, 763, 124, 489, 371, 313, 273, 753, 815, 879, 573, 629, 508, 624, 350, 89, 936, 522, 673, 889, 275, 757, 349, 412, 343, 287, 897, 129, 190, 321, 391, 611, 240, 595, 332, 155, 761, 6, 95, 743, 342, 176, 107, 262, 245, 111, 811, 896, 232, 182, 31, 843, 891, 78, 189, 964, 135, 479, 14, 272, 512, 590, 24, 56, 52, 484, 988, 365, 555, 198, 76, 192, 238, 511, 320, 94, 228, 603, 963, 794, 75, 237, 337, 894, 137, 601, 143, 164, 154, 440, 50, 145, 420, 123, 941, 535, 842, 432, 790, 294, 664, 375, 764, 121, 333, 475, 252, 860, 372, 808, 357, 533, 909, 773, 53, 46, 326, 26, 680, 569, 495, 892, 861, 393, 351, 330, 870, 858, 907, 235, 709, 902, 827, 816, 118, 100, 345, 961, 698, 562, 139, 615, 243, 924, 149, 474, 43, 583, 230, 962, 610, 987, 574, 445, 748, 425, 331, 379, 225, 166, 253, 102, 358, 755, 684, 162, 836, 431, 422, 413, 276, 552, 157, 396, 654, 284, 110, 838, 599, 471, 112, 598, 933, 397, 490, 42, 521, 934, 659, 513, 912, 967, 687, 818, 768, 209, 766, 202, 730, 11, 636, 410, 19, 619, 150, 310, 553, 325, 27, 504, 7, 250, 120, 671, 444, 847, 529, 233, 20, 8, 968, 523, 826, 778, 951, 985, 839, 387, 950, 442, 295, 200, 977, 196, 309, 884, 70, 644, 706, 797, 947, 403, 596, 1, 848, 908, 645, 406, 738, 459, 903, 817, 825, 953, 96, 530, 64, 661, 732, 593, 678, 681, 656, 423, 646, 675, 151, 604, 981, 323, 801, 719, 361, 898, 125, 340, 700, 199, 226, 482, 424, 187, 288, 353, 524, 328, 103, 837, 954, 608, 404, 980, 158, 887, 483, 734, 221, 647, 587, 501, 558, 268, 798, 214, 557, 428, 650, 995, 473, 708, 942, 376, 307, 346, 152, 146, 85, 217, 754, 500, 747, 49, 921, 575, 690, 465, 791, 492, 701, 640, 869, 637, 395, 772, 740, 218, 991, 758, 641, 622, 260, 29, 115, 972, 862, 463, 34, 746, 813, 971, 44, 377, 831, 167, 904, 175, 551, 304, 616, 383, 464, 458, 73, 949, 715, 293, 140, 767, 438, 733, 286, 194, 67, 77, 266, 485, 899, 945, 745, 244, 363, 854, 565, 84, 721, 299, 620, 335, 472, 303, 223, 443, 613, 792, 138, 682, 279, 457, 657, 705, 261, 914, 347, 990, 329, 41, 285, 82, 45, 528, 806, 605, 957, 779, 660, 486, 973, 517, 696, 212, 919, 916, 602, 821, 434, 263, 850, 239, 756, 634, 824, 241, 356, 516, 3, 122, 165, 0, 514, 509, 814, 714, 222, 362, 271, 863, 86, 9, 254, 116, 662, 974, 744, 318, 699, 249, 302, 793, 359, 105, 36, 789, 929, 688, 153, 752, 136, 691, 539, 695, 906, 355, 567, 314, 655, 178, 255, 206, 623, 651, 297, 742, 588, 494, 983, 81, 751, 875, 40, 643, 267, 88, 213, 844, 538, 283, 829, 849, 334, 210, 630, 367, 370, 866, 246, 857, 716, 771, 802, 305, 452, 581, 134, 127, 205, 881, 292, 462, 352, 668, 117, 30, 227, 880, 156, 37, 717, 290, 952, 589, 269, 937, 585, 799, 13, 710, 626, 917, 378, 505, 614, 98, 503, 208, 976, 712, 873, 648, 448, 257, 17, 63, 360, 895, 550, 686, 635, 515, 617, 437, 874, 609, 702, 204, 429};
    }

    // for (size_t i = 0; i < cap; i++)
    // {
    //     for (size_t j = 0; j < cap; j++)
    //     {
    //         fprintf(stderr,"%zu,%zu,%.2f\n",i,j,calcDistance(seqs[foo[i]], seqs[foo[j]]));
    //     }
    // }
    size_t clus = tree.first_insert(seqs[foo[0]], insertionList, foo[0]);
    for (size_t i = 1; i < cap; i++)
    {
        // fprintf(stderr, "inserting %zu\n", foo[i]);
        size_t clus = tree.insert(seqs[foo[i]], insertionList, foo[i]);
        // clusters[foo[i]] = tree.findAncestor(clus);
        clusters[foo[i]] = clus;
    }
    // for (size_t i = 0; i < cap; i++)
    // {
    //     clusters[foo[i]] = tree.superCluster(clusters[foo[i]]);
    // }
    fprintf(stderr, "\n\n\n\n");
    tree.printTreeJson(stderr);

    // size_t n = 102;
    // double threshold = 1.5;
    // tree.forceSplit(n, insertionList, threshold);

    // for (size_t child : tree.childLinks[n])
    // {
    //     if (tree.isBranchNode[child] && tree.priority[child] > threshold)
    //     {
    //         tree.forceSplit(child, insertionList, threshold);
    //     }
    // }

    // size_t best = tree.findNearest(seqs[15], 19, 0);
    // fprintf(stderr,"***%zu\n",best);

    // fprintf(stderr, "%f, %f \n", tree.means[48].first, tree.means[48].second);

    for (size_t run = 0; run < iteration; run++)
    {
        tree.printTreeJson(stdout);
        fprintf(stderr, "Iteration %zu\n", run);
        tree.prepReinsert();
        for (size_t i = 0; i < cap; i++)
        {
            size_t clus = tree.reinsert(seqs[foo[i]], foo[i]);

            // fprintf(stderr, "\n found %zu at %zu\n", foo[i], clus);
            // clusters[foo[i]] = tree.findAncestor(clus);
            clusters[foo[i]] = clus;
        }

        tree.trim(insertionList[insertionList.size() - 1] - 1);
    }

    FILE *pFile;
    if (random_)
    {
        pFile = fopen("tree-test-r.csv", "w");
    }
    else
    {
        pFile = fopen("tree-test.csv", "w");
    }
    fprintf(pFile, "node,parent,level,isBranch,priority,x,y\n");
    tree.printSubTreeMatrices(pFile, seqs, 0);

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

    // // Recursively destroy all locks
    // tree.destroyLocks();

    tree.printTreeJson(stdout);
    tree.destroyTree();
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
    vector<pair<float, float>> centroids = generateCentroids(clusCount, 5);

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
    split_threshold = 5;
    stay_threshold = 1;
    minimiser_match_threshold = 4;
    partree_capacity = 10000;

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
    split_node_threshold = split_threshold / 2;
    // split_node_threshold = stay_threshold + 2;

    fprintf(stderr, "split threshold: %.2f\n", split_threshold);
    fprintf(stderr, "stay threshold: %.2f\n", stay_threshold);
    fprintf(stderr, "split_node_threshold threshold: %.2f\n", split_node_threshold);

    if (args.minimiser_match_given)
    {
        minimiser_match_threshold = args.minimiser_match_arg;
        fprintf(stderr, "minimiser_match threshold: %zu\n", minimiser_match_threshold);
    }

    iteration = args.iteration_arg;

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
