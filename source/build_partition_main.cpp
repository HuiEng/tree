#include <seqan3/search/views/partition_multi_hash.hpp>
#include <seqan3/search/views/partition_multi.hpp>
#include "minimiser.hpp"
#include "bloom_filter.hpp"
// #include "part_ktree.hpp"
#include <filesystem>
#include "build_partition_main_cmdline.hpp"

using namespace std;
namespace fs = std::filesystem;

static build_partition_main_cmdline args; // Command line switches and arguments
static uint8_t kmerLength = 9;            // Kmer length
static uint32_t windowLength = 50;        // window length
static uint32_t step_size = windowLength; // window length
size_t minimiser_size = 3;
size_t bf_element_cnt = 1000;
bool debug = false;
bool compressReads = false;
bool multipleOut = false;
string outfile = "";
size_t chunkRatio = 1;

// void writeInt(std::ostream &os, unsigned long long int i)
// {
//     os.write(reinterpret_cast<const char *>(&i), sizeof(i));
// }

template <typename view>
void getMinimisers(view minimiser_view, bloom_parameters parameters, string filename, ofstream &wf)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    size_t max = 8;
    size_t end = -1;

    bloom_filter bf(parameters);

    if (compressReads)
    {
        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
            }
        }
        bf.print(wf);
        bf.clear();
    }
    else
    {
        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
            }
            bf.print(wf);
            bf.clear();
        }
    }
}

template <typename view>
void getPartitionMinimisers(view minimiser_view, bloom_parameters parameters, string filename, ofstream &wf)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    // ofstream outfile(outname);
    size_t max = 8;
    // ofstream wf(outname, ios::out | ios::binary);
    size_t end = -1;

    bloom_filter bf(parameters);
    // writeInt(wf, bf.table_size());

    if (debug)
    {
        size_t i = 0;

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            // fprintf(stdout, ">\n");
            cout << ">0|" << i << "|0\n";
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                    cout << hashToMer_str(kmerLength, hash) << ";";
                }
                cout << ",";
                bf.print(wf);
                // bf.printBFIdx(stderr);

                // ///// debug
                // bf.print();

                bf.clear();
            }
            cout << "\n";

            // fprintf(stderr,">s\n");

            // end of seq flag, print empty bf
            bf.print(wf);
            i++;
        }
    }
    else if (compressReads)
    {
        auto first_seq = (*file_in.begin()).sequence();
        auto windows = first_seq | minimiser_view;
        vector<bloom_filter> wBFL;
        for (auto &&hashes : windows)
        {
            bloom_filter temp(parameters);
            wBFL.push_back(temp);
        }

        for (auto &[seq, id, qual] : file_in)
        {
            size_t i = 0;
            // fprintf(stdout, ">\n");
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    wBFL[i].insert(hash);
                }
                i++;
            }
        }

        // print the list
        for (bloom_filter wbf : wBFL)
        {
            wbf.print(wf);
        }
        // end of seq flag, print empty bf
        bf.print(wf);
    }

    else
    {

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            auto result = seq  | std::views::reverse | seqan3::views::complement | minimiser_view;
            auto it = result.begin();

            // fprintf(stdout, ">\n");
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
                for (size_t hash : *it)
                {
                    bf.insert(hash);
                }
                it++;
                bf.print(wf);
                // bf.printBFIdx(stderr);

                // ///// debug
                // bf.print();

                bf.clear();
            }

            // fprintf(stderr,">s\n");

            // end of seq flag, print empty bf
            bf.print(wf);
        }
    }
    // wf.close();
}



template <typename view>
void compressPartitionMinimisers(view minimiser_view, bloom_parameters parameters, string filename, ofstream &wf)
{
    seqan3::sequence_file_input<dna4_traits> file_in{filename};
    // ofstream outfile(outname);
    size_t max = 8;
    // ofstream wf(outname, ios::out | ios::binary);
    size_t end = -1;

    bloom_filter bf(parameters);
    // writeInt(wf, bf.table_size());

    if (debug)
    {
        size_t i = 0;

        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            // fprintf(stdout, ">\n");
            cout << ">0|" << i << "|0\n";
            size_t n = 0;
            auto result = seq  | std::views::reverse | seqan3::views::complement | minimiser_view ;
            // cout << "***"<<std::distance(result.begin(), result.end()) <<"\n";
             cout << "***\n";
             auto it = result.begin();
            // seqan3::debug_stream << *it << '\n';
            // seqan3::debug_stream << *(result.end()) << '\n';

            // for (size_t i = 0;i<hashes.size();i++)
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                    cout << hashToMer_str(kmerLength, hash) << ";";
                }
                cout << "\n";
                for (size_t hash : *it)
                {
                    bf.insert(hash);
                    cout << hashToMer_str(kmerLength, hash) << ";";
                }
                cout << "\n";
                seqan3::debug_stream << hashes << '\n';
                seqan3::debug_stream << *it << '\n';
                it++;
                n++;
                if (n==chunkRatio){
                    cout << "###########\n";
                    bf.print(wf);
                    bf.clear();
                    n=0;
                }
            }
            if (n!=0){
                cout << "@@@\n";
                bf.print(wf);
                bf.clear();
                n=0;
            }
            cout << "\n";
            // end of seq flag, print empty bf
            bf.print(wf);
            i++;
        }
    }
    else
    {
        size_t n = 0;
        // Retrieve the sequences and ids.
        for (auto &[seq, id, qual] : file_in)
        {
            auto result = seq  | std::views::reverse | seqan3::views::complement | minimiser_view;
            auto it = result.begin();
            // fprintf(stdout, ">\n");
            for (auto &&hashes : seq | minimiser_view)
            {
                for (size_t hash : hashes)
                {
                    bf.insert(hash);
                }
                for (size_t hash : *it)
                {
                    bf.insert(hash);
                }
                it++;
                n++;
                if (n==chunkRatio){
                    bf.print(wf);
                    bf.clear();
                    n=0;
                }
            }

            if (n!=0){
                bf.print(wf);
                bf.clear();
                n=0;
            }    

            // end of seq flag, print empty bf
            bf.print(wf);
        }
    }
    // wf.close();
}


void doWork(ofstream &wf, bloom_parameters parameters, string inputFile)
{
    // ofstream wf(outfile, ios::out | ios::binary);
    // bloom_filter bf(parameters);
    // writeInt(wf, bf.table_size());

    if (args.canonical_arg)
    {
        auto partition_view = seqan3::views::partition_multi_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
                                                                  seqan3::window_size{windowLength},
                                                                  seqan3::minimiser_size{minimiser_size},
                                                                  seqan3::step_size{step_size},
                                                                  seqan3::seed{0});
        getPartitionMinimisers(partition_view, parameters, inputFile, wf);
    }
    else
    {
        // to get minimisers with w=8,k=4
        // input param for the minimiser view is calculated by: window size - k-mer size + 1, here: 8 - 4 + 1 = 5)
        size_t temp = windowLength - kmerLength + 1;
        auto partition_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::partition_multi(temp, kmerLength, minimiser_size, step_size);
// auto partition_view = seqan3::views::partition_multi_hash(seqan3::shape{seqan3::ungapped{kmerLength}},
//                                                                   seqan3::window_size{windowLength},
//                                                                   seqan3::minimiser_size{minimiser_size},
//                                                                   seqan3::step_size{step_size},
//                                                                   seqan3::seed{0});
        if (args.toSingle_arg)
        {
            getMinimisers(partition_view, parameters, inputFile, wf);
        }
        else if (chunkRatio != 1)
        {
            // double chunkRatio = windowLength*1.0/step_size;
            // windowLength = step_size;
            // size_t temp = windowLength - kmerLength + 1;
            // auto partition_view = seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{kmerLength}}) | seqan3::views::partition_multi(temp, kmerLength, minimiser_size, step_size);

            compressPartitionMinimisers(partition_view, parameters, inputFile, wf);
        }
        else
        {
            getPartitionMinimisers(partition_view, parameters, inputFile, wf);
        }
    }

    // wf.close();
}


int build_partition_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    //

    debug = args.debug;

    if (args.kmer_given)
        kmerLength = args.kmer_arg;
    if (args.window_given)
        windowLength = args.window_arg;

    if (kmerLength > windowLength)
    {
        fprintf(stderr, "Error: kmer length must be smaller or equal to window length\n");
        return 1;
    }

    if (args.size_arg > windowLength - kmerLength + 1)
    {
        fprintf(stderr, "Error: number of minimisers per window must be smaller than %d (w - k + 1)\n", windowLength - kmerLength + 1);
        return 1;
    }

    if (args.element_given)
        bf_element_cnt = args.element_arg;

    bloom_parameters parameters;
    // How many elements roughly do we expect to insert?
    parameters.projected_element_count = bf_element_cnt;

    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = 0.001; // 1 in 10000

    // Simple randomizer (optional)
    parameters.random_seed = 0xA5A5A5A5;
    parameters.maximum_number_of_hashes = 1;

    if (!parameters)
    {
        std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
        return 1;
    }
    parameters.compute_optimal_parameters();

    compressReads = args.compress_arg;
    if (compressReads)
    {
        std::cout << "Compressing Reads" << std::endl;
    }
    multipleOut = args.multiple_arg;
    if (multipleOut)
    {
        std::cout << "multipleOut" << std::endl;
    }

    if (args.size_given)
    {
        minimiser_size = args.size_arg;
    }
    fprintf(stderr, "Partition - Generating %zu minimisers per window...\n", minimiser_size);
    fprintf(stderr, "kmerLength= %u, windowLength = %u\n", kmerLength, windowLength);

    string inputFile = args.input_arg;
    size_t firstindex = inputFile.find_last_of("/") + 1;
    size_t lastindex = inputFile.find_last_of(".");
    outfile = inputFile.substr(firstindex, lastindex - firstindex);
    string buffer = "";
    char bufferArr[50];
    if (args.element_given)
    {
        bf_element_cnt = args.element_arg;
    }
    
    if (args.step_given)
    {
        step_size = args.step_arg;
        sprintf(bufferArr, "-k%u-w%u-s%zu-b%zu--step%u", kmerLength, windowLength, minimiser_size, bf_element_cnt, step_size);
    }
    else
    {
        step_size = windowLength;
        sprintf(bufferArr, "-k%u-w%u-s%zu-b%zu", kmerLength, windowLength, minimiser_size, bf_element_cnt);
    }
    buffer = buffer + bufferArr;
    if (args.toSingle_arg)
    {
        buffer = buffer + "-single";
    }

    double temp = windowLength*1.0/step_size;
    chunkRatio = temp;
    windowLength = step_size;
    fprintf(stderr,"chunkRatio %zu\n", chunkRatio);

    if (args.folder_arg)
    {
        if (args.output_given)
        {
            outfile = args.output_arg;
        }
        else
        {
            fprintf(stderr, "Please provide output file path\n");
            return 0;
        }
        
        outfile = outfile + buffer;
        bloom_filter bf(parameters);
        outfile = outfile + ".bin";
        ofstream wf(outfile, ios::out | ios::binary);
        writeInt(wf, bf.table_size());

        string line;
        size_t cnt = 0;

        // Read from the text file
        ifstream listStream(inputFile.c_str());
        while (getline(listStream, line))
        {
            doWork(wf, parameters, line);
            cnt++;
        }
        listStream.close();
        wf.close();
        fprintf(stderr,"Processed %zu files, signatureSize %zu\n", cnt, bf.table_size());
        fprintf(stderr,"lastFile to %s signatureSize %zu\n", line.c_str(), bf.table_size());
        fprintf(stderr,"output to %s \n", outfile.c_str());
        return 0;

        // string delimiter = "/*";
        // string folder = inputFile.substr(0, inputFile.find(delimiter));
        // string ext = inputFile.substr(inputFile.find(delimiter) + delimiter.size(), inputFile.size() - 1);
        // fprintf(stderr, "Reading folder %s\n", folder.c_str());

        // outfile = outfile + buffer;

        // bloom_filter bf(parameters);

        // if (multipleOut)
        // {
        //     size_t i = 0;
        //     for (const auto &entry : fs::directory_iterator(folder))
        //     {
        //         if (entry.path().extension() == ext)
        //         {
        //             ofstream wf(outfile + "_" + to_string(i) + ".bin", ios::out | ios::binary);
        //             writeInt(wf, bf.table_size());
        //             doWork(wf, parameters, entry.path());
        //             wf.close();
        //             i++;
        //         }
        //     }
        // }
        // else
        // {
        //     outfile = outfile + ".bin";
        //     ofstream wf(outfile, ios::out | ios::binary);
        //     writeInt(wf, bf.table_size());
        //     for (const auto &entry : fs::directory_iterator(folder))
        //     {
        //         if (entry.path().extension() == ext)
        //         {
        //             // cout << entry.path().stem().string() << '\n';
        //             doWork(wf, parameters, entry.path());
        //         }
        //     }
        //     wf.close();
        //     return 0;
        // }
    }

    else
    {
        if (args.output_given)
        {
            outfile = args.output_arg;
        }
        outfile = outfile + buffer + ".bin";
        ofstream wf(outfile, ios::out | ios::binary);
        bloom_filter bf(parameters);
        writeInt(wf, bf.table_size());
        doWork(wf, parameters, inputFile);
        wf.close();
    }
    return 0;
}
