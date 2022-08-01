#include <random>
#include "tree_main_cmdline.hpp"
#include "self_tree.hpp"
typedef self_tree tree_type;

using namespace std;

static tree_main_cmdline args; // Command line switches and arguments

void outputClusters(FILE *pFile, const vector<tuple<size_t, size_t>> &clusters)
{
    fprintf(pFile, "seqID,cluster,ancestor\n");
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu,%llu\n",
                static_cast<unsigned long long>(sig),
                static_cast<unsigned long long>(get<0>(clusters[sig])),
                static_cast<unsigned long long>(get<1>(clusters[sig])));
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

vector<tuple<size_t, size_t>> clusterSignatures(const vector<vector<vector<cell_type>>> &seqs)
{
    size_t seqCount = seqs.size();
    // seqCount = 26;
    vector<tuple<size_t, size_t>> clusters(seqCount);
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
    // obtain a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // shuffle(foo.begin(), foo.end(), std::default_random_engine(seed));
    // vector<size_t> foo{212,677,185,710,785,1035,512,29,690,594,929,1016,801,683,1065,211,390,789,945,408,806,687,302,736,66,260,649,507,771,100,72,417,599,1010,981,439,563,484,291,969,11,144,328,362,530,588,906,403,540,146,832,96,351,149,319,254,104,986,971,777,944,158,639,414,421,189,889,611,52,863,276,161,850,941,853,22,962,565,486,531,591,124,127,210,247,1036,815,516,749,514,1089,462,696,31,511,536,650,763,398,379,520,1019,356,823,373,1063,384,267,557,350,263,383,664,891,141,358,446,692,852,866,36,98,34,524,570,139,103,528,663,503,290,933,521,1072,975,1077,74,829,908,265,502,497,1015,54,874,917,85,821,416,655,733,376,91,345,722,883,122,916,743,938,270,946,605,830,313,134,130,812,1062,756,553,1004,231,355,1055,729,459,487,781,311,208,926,495,921,934,930,1030,651,1040,839,457,537,6,509,741,389,949,136,714,387,1027,489,798,175,219,242,256,545,898,413,152,207,773,300,449,407,738,108,296,593,264,21,622,162,62,887,814,105,79,787,166,816,470,1061,835,499,236,970,222,221,595,151,586,783,884,577,834,923,20,28,619,43,548,755,1042,450,478,485,227,55,75,615,238,112,993,107,713,431,653,673,168,602,1082,286,732,872,145,435,253,643,1031,361,1068,1014,468,542,626,498,354,1029,634,332,697,871,120,997,967,335,724,681,609,966,1079,440,183,106,554,578,434,566,510,747,337,341,752,474,123,942,841,186,372,70,876,125,766,39,129,192,2,775,802,523,323,196,334,1025,411,438,625,49,339,700,527,885,19,481,188,796,1086,117,652,38,840,367,1000,948,243,1081,466,678,228,460,657,564,301,271,342,92,712,76,456,644,241,904,86,237,1032,371,689,820,453,613,282,568,277,953,892,534,1003,479,894,709,142,828,465,734,215,95,306,826,443,27,206,680,1083,232,508,73,786,1091,119,138,719,285,1070,377,620,10,193,7,546,1018,702,667,58,541,701,89,327,817,629,261,784,116,698,896,688,860,101,396,668,441,574,225,859,430,169,632,601,424,1041,1007,288,148,822,472,454,137,980,529,659,195,268,297,899,312,581,4,903,60,999,1075,458,121,153,233,589,473,746,910,165,316,513,394,381,776,645,173,1071,562,1005,375,71,932,555,308,150,772,155,87,401,630,331,901,349,99,199,24,272,902,410,273,329,1006,984,654,17,220,287,600,30,596,333,269,352,661,726,936,496,190,637,476,842,167,1090,675,404,32,1046,5,847,492,110,867,406,543,811,631,1023,201,244,303,3,275,1059,1066,585,788,880,77,990,234,205,760,735,204,235,131,1020,405,873,539,592,448,567,556,289,368,1060,51,480,587,606,935,579,769,656,610,983,365,742,177,658,865,704,624,69,895,452,477,744,59,877,132,1056,957,279,616,44,344,793,517,538,836,685,309,745,778,346,1052,897,451,395,1050,703,1044,988,1034,156,638,409,40,179,33,283,751,397,1087,53,184,708,194,780,42,315,768,900,45,322,647,571,419,386,790,1053,392,956,947,359,608,711,65,258,172,176,181,582,266,964,753,483,612,378,532,280,257,500,83,813,764,730,1008,918,230,870,809,559,939,494,427,444,1064,490,757,914,535,0,202,774,318,958,799,111,731,491,133,135,976,363,607,353,203,695,209,180,164,1073,1080,366,475,447,295,940,875,795,41,1048,982,991,170,597,576,636,862,1021,779,428,580,1,37,791,23,825,429,627,728,996,218,293,623,81,420,864,274,515,94,1076,759,721,461,869,82,330,455,157,590,909,522,126,223,837,748,699,937,550,357,573,818,118,423,26,886,84,12,109,720,762,868,18,213,754,171,340,525,679,326,248,469,348,977,959,1013,694,603,229,1085,992,1017,952,1038,851,246,418,501,214,46,400,114,432,666,464,385,973,640,262,672,575,633,154,831,1037,890,968,642,504,370,216,97,393,974,665,78,187,737,670,245,493,1047,648,1002,9,569,907,669,995,422,950,794,255,979,128,888,433,425,294,922,198,739,314,321,347,217,572,805,765,88,800,912,90,682,163,920,707,716,833,547,845,723,61,1069,505,924,558,115,965,93,1074,718,925,767,843,374,1088,64,854,304,761,552,382,770,249,808,369,25,943,662,671,415,855,913,343,758,1001,824,67,1009,445,684,693,931,1024,676,8,725,893,985,13,963,191,810,846,583,519,1051,467,686,621,844,560,226,292,338,1045,250,298,336,614,113,827,544,402,803,1043,972,928,882,1057,807,200,15,240,961,388,1026,1084,750,955,436,50,182,437,838,1012,1058,617,717,533,259,251,706,646,1067,463,878,674,989,364,252,715,63,978,960,281,391,551,56,160,147,797,14,911,488,518,857,740,380,915,305,1078,905,80,16,35,927,482,1011,628,224,804,299,954,727,399,604,549,998,819,284,881,660,1039,239,951,506,278,140,197,691,174,68,47,360,858,994,1022,471,1049,442,48,598,792,584,412,310,861,526,102,324,143,856,919,307,705,426,178,879,1054,320,1028,635,782,618,1033,561,57,317,159,849,325,848,987};

    // Insert first 1 nodes single-threaded
    for (size_t i = 0; i < firstNodes; i++)
    {
        size_t clus = tree.first_insert(seqs[foo[i]], insertionList, foo[i]);
        clusters[foo[i]] = make_tuple(clus, clus);
        // fprintf(stderr, ">>%zu\n", foo[i]);
        // dbgPrintSignatureIdx(stderr, seqs[foo[i]]);
    }

    for (size_t i = firstNodes; i < seqCount; i++)
    {
        // fprintf(stdout, "inserting %zu", i);
        size_t clus = tree.insert(seqs[foo[i]], insertionList, foo[i]);
        // clusters[foo[i]] = tree.findAncestor(clus);
        clusters[foo[i]] = make_tuple(clus, tree.findAncestor(clus));
        fprintf(stderr, ">>%zu\n", foo[i]);
        // dbgPrintSignatureIdx(stderr, seqs[foo[i]]);
    }

    // tree.removeSingleton(clusters, insertionList);

    // //search
    // for (size_t i = 0; i < seqCount; i++)
    // {
    //     // fprintf(stdout, "inserting %zu", foo[i]);
    //     size_t clus = tree.search(seqs[foo[i]], foo[i]);
    //     // clusters[foo[i]] = tree.findAncestor(clus);
    //     clusters[foo[i]] = make_tuple(clus, tree.findAncestor(clus));
    // }

    // Recursively destroy all locks
    tree.destroyLocks();

    tree.printTreeJson(stdout);

    return clusters;
}

int tree_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    //?
    if (!args.input_given)
    {
        cout << "No input and/or query given! Exiting...\n";
        return 0;
    }

    if (args.split_threshold_given)
    {
        split_threshold = args.split_threshold_arg;
        fprintf(stderr, "split threshold: %.2f\n", split_threshold);
    }

    if (args.stay_threshold_given)
    {
        stay_threshold = args.stay_threshold_arg;
        fprintf(stderr, "stay threshold: %.2f\n", stay_threshold);
    }

    if (args.minimiser_match_given)
    {
        minimiser_match_threshold = args.minimiser_match_arg;
        fprintf(stderr, "minimiser_match threshold: %zu\n", minimiser_match_threshold);
    }

    string inputFile = args.input_arg;

    vector<vector<vector<cell_type>>> seqs = readPartitionBF(inputFile);
    fprintf(stderr, "Loaded signatures...\n");

    signatureSize = seqs[0][0].size();
    fprintf(stderr, "Building Signature...\n");
    default_random_engine rng;
    vector<tuple<size_t, size_t>> clusters = clusterSignatures(seqs);

    fprintf(stderr, "writing output...\n");

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    auto fileName = rawname + "-s" + to_string((int)(stay_threshold * 100)) + "-l" + to_string((int)(split_threshold * 100));
    FILE *pFile = fopen((fileName + ".txt").c_str(), "w");
    outputClusters(pFile, clusters);

    return 0;
}
