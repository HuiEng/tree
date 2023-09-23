// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "similarity_main_cmdline.hpp"
#include "similarity.hpp"
static similarity_main_cmdline args; // Command line switches and arguments

void allKmersBatch(FILE *pFile, const string inputFile)
{
    // pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
    size_t idx = 0;
    vector<cell_type> seqs;
    signatureSize = readSignatures(inputFile, seqs, batch, idx);

    fprintf(pFile, "i,j,similarity\n");
    calcAllSimilarityKmers(pFile, seqs);

    size_t start = idx;
    vector<cell_type> temp;
    vector<cell_type> seqsA;
    size_t startA = start;
    size_t idxA = 0;
    bool changed = false;

    do
    {
        readSignatures(inputFile, temp, batch, idx);
        size_t offset = (idx - temp.size()) / signatureSize;
        calcAllSimilarityKmers(pFile, temp, offset);
        calcAllSimilarityKmersBatch(pFile, seqs, temp, idxA / signatureSize, offset);

        if (!changed)
        {
            seqsA = temp;
            startA = start;
            changed = true;
        }
        start = idx;

    } while (temp.size() > 0);

    size_t seqCount = idx;
    while (seqsA.size() > 1)
    {
        seqs = seqsA;
        idx = startA + seqsA.size();
        start = idx;
        idxA = startA;

        changed = false;

        do
        {
            readSignatures(inputFile, temp, batch, idx);
            size_t offset = (idx - temp.size()) / signatureSize;
            calcAllSimilarityKmersBatch(pFile, seqs, temp, idxA / signatureSize, offset);

            if (!changed)
            {
                seqsA = temp;
                startA = start;
                // idxA = idx;
                changed = true;
            }
            start = idx;

        } while (temp.size() > 0);
    }
    seqCount = seqCount / signatureSize;
    fprintf(stderr, "Loaded %zu seqs..\n", seqCount);
}

void outputClusters(FILE *pFile, const vector<size_t> &clusters)
{
    for (size_t sig = 0; sig < clusters.size(); sig++)
    {
        fprintf(pFile, "%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
    }
}

vector<vector<size_t>> compressClusterList(const string file)
{
    string line;
    ifstream clusterStream(file);

    size_t id = 0;
    size_t clus = 0;
    vector<size_t> clusters;

    while (getline(clusterStream, line))
    {
        sscanf(line.c_str(), "%zu,%zu", &id, &clus);
        clusters.push_back(clus);
    }
    clusterStream.close();

    unordered_map<size_t, size_t> remap;
    vector<vector<size_t>> clusterGroup;
    size_t i = 0;
    for (size_t &clus : clusters)
    {
        if (remap.count(clus))
        {
            clus = remap[clus];
            clusterGroup[clus].push_back(i);
        }
        else
        {
            size_t newClus = remap.size();
            remap[clus] = newClus;
            clus = newClus;
            clusterGroup.push_back({i});
        }
        i++;
    }
    fprintf(stderr, "Output %zu clusters\n", remap.size());
    FILE *cFile = fopen("test.txt", "w");
    outputClusters(cFile, clusters);
    return clusterGroup;
}

template <typename T>
void insertVecRange(vector<T> &to, vector<T> &from)
{
    to.insert(to.end(), from.begin(), from.end());
}

template <typename funcType, typename sigType>
void sepCluster(const string file,
                funcType simFunc,
                vector<sigType> &seqs, size_t mul = 1)
{
    string rawname = "test";

    vector<vector<size_t>> clusterGroup = compressClusterList(file);
    size_t offset = 0;
    size_t offset_old = 0;
    size_t clu = 0;
    for (vector<size_t> cluster : clusterGroup)
    {
        FILE *pFile = fopen(("intra" + to_string(clu) + "-all_sim.txt").c_str(), "w");
        fprintf(pFile, "i,j,similarity\n");
        vector<sigType> temp;
        for (size_t id : cluster)
        {
            if (mul == 1)
            {
                temp.push_back(seqs[id]);
            }
            else
            {
                sigType *signature = &seqs[id * mul];
                temp.insert(temp.end(), signature, signature + signatureSize);
            }
            offset++;
        }
        fprintf(stderr, "Loaded %zu seqs...in cluster %zu\n", temp.size() / mul, clu);
        clu++;
        simFunc(pFile, temp, offset_old);
        offset_old = offset;
    }
}

int similarity_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string inputFile = args.bf_input_arg;

    if (args.threshold_given)
        minimiser_match_threshold = args.threshold_arg;

    if (args.skip_arg)
    {
        skip = true;
    }

    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    string rawname = inputFile.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }
    FILE *pFile;

    if (args.cluster_given)
    {
        if (args.all_kmer_arg)
        {
            vector<cell_type> seqs;
            signatureSize = readSignatures(inputFile, seqs);
            fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
            sepCluster(args.cluster_arg, calcAllSimilarityKmers, seqs, signatureSize);
        }
        else
        {
            vector<seq_type> seqs = readPartitionBF(inputFile, signatureSize);
            fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());
            sepCluster(args.cluster_arg, &calcAllSimilarityWindow, seqs);
        }
    }
    else if (args.batch_given)
    {
        skip = true;
        batch = args.batch_arg;
        if (args.all_kmer_arg)
        {
            pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
            allKmersBatch(pFile, inputFile);
            // vector<vector<cell_type>> seqs_batch = readSignaturesBatch(inputFile, batch, signatureSize);
            // size_t temp = seqs_batch.size() - 1;
            // size_t seqCount = temp * batch + seqs_batch[temp].size() / signatureSize;
            // fprintf(stderr, "Loaded %zu seqs...\n", seqCount);
            // batchSimFunction(pFile, seqs_batch, &calcAllSimilarityKmers, &calcAllSimilarityKmersBatch);
        }
        else
        {
            // vector<vector<vector<vector<cell_type>>>> seqs_batch = readPartitionBFBatch(inputFile, batch, signatureSize);
            // size_t temp = seqs_batch.size() - 1;
            // size_t seqCount = temp * batch + seqs_batch[temp].size();
            // fprintf(stderr, "Loaded %zu seqs...\n", seqCount);

            if (args.local_arg)
            {
                pFile = fopen((rawname + "local_sim.txt").c_str(), "w");
                batching(inputFile, pFile, &calcAllSimilarityLocal, &calcAllSimilarityBatch);
            }
            else if (args.global_arg)
            {
                pFile = fopen((rawname + "-global_sim.txt").c_str(), "w");
                batching(inputFile, pFile, &calcAllSimilarityGlobal, &calcAllSimilarityGlobalBatch);
                // batchSimFunction(pFile, seqs_batch, &calcAllSimilarityGlobal, &calcAllSimilarityGlobalBatch);
            }
            else
            {
                char buffer[50];
                sprintf(buffer, "-t%zu-window_sim.txt", minimiser_match_threshold);
                pFile = fopen((rawname + buffer).c_str(), "w");
                batching(inputFile, pFile, &calcAllSimilarityWindow, &calcAllSimilarityLocalBatch);
            }
        }
    }
    else
    {
        if (args.all_kmer_arg)
        {
            pFile = fopen((rawname + "-all_sim.txt").c_str(), "w");
            vector<cell_type> seqs;
            if (args.multiple_arg)
            {
                signatureSize = readList(inputFile, seqs);
            }
            else
            {
                signatureSize = readSignatures(inputFile, seqs);
            }

            fprintf(stderr, "Loaded %zu seqs...\n", seqs.size() / signatureSize);
            fprintf(pFile, "i,j,similarity\n");

            calcAllSimilarityKmers(pFile, seqs);
        }
        else
        {
            vector<seq_type> seqs;
            if (args.multiple_arg)
            {
                seqs = readListPartitionBF(inputFile, signatureSize);
            }
            else
            {
                seqs = readPartitionBF(inputFile, signatureSize);
            }

            if (signatureSize == 0)
            {
                fprintf(stderr, "Something is wrong with the input data, please generate signature with diff params\n");
                return 0;
            }

            // fprintf(stderr,"done\n" );
            fprintf(stderr, "Loaded %zu seqs...\n", seqs.size());

            // calcAllSetBits(sigs);
            if (args.local_arg)
            {
                pFile = fopen((rawname + "-local_sim.txt").c_str(), "w");
                fprintf(pFile, "i,j,similarity\n");

                calcAllSimilarityLocal(pFile, seqs);
            }
            else if (args.global_arg)
            {
                pFile = fopen((rawname + "-global_sim.txt").c_str(), "w");
                fprintf(pFile, "i,j,similarity\n");

                calcAllSimilarityGlobal(pFile, seqs);
            }
            else
            {

                char buffer[50];
                sprintf(buffer, "-t%zu-window_sim.txt", minimiser_match_threshold);
                pFile = fopen((rawname + buffer).c_str(), "w");
                fprintf(pFile, "i,j,similarity\n");
                calcAllSimilarityWindow(pFile, seqs);
            }
        }
    }

    // for (int i = 0; i < sigs.size(); i += signatureSize)
    // {
    //     printBF(&sigs[i]);
    // }

    // for (auto sig:sigs){
    //     fprintf(stdout,"%zu\n",sig);
    // }

    return 0;
}
