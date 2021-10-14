// Water_order.cpp : This file contains the 'main' function. Program execution begins and ends there.
/* Oriental tetrahedral order parameter calculator
 Written by Shoubhik Raj Maiti, Aug 2021 */

/* Reading of chemical structure files is done with Chemfiles by Guillaume Fraux
 Command line argument parsing is done with TCLAP by Mike Smoot and Daniel Aarno */

// 

#include <chemfiles.hpp>
#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <array>
#include <tclap/CmdLine.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include <cfenv>
#include <voro++.hh>

// when using MSVC++, include basetsd.h for SSIZE_T because it's ancient OpenMP 2.0
// does not support unsigned index for OpenMP for loops. Intel C++ also defines
//  _MSC_VER on Windows but they support OpenMP > 4.0 so check for them explicitly.
// In those compilers, there is no need to use unsigned index!!
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
 #include <basetsd.h>
 typedef SSIZE_T openmp_index_type;
#else
 typedef size_t openmp_index_type;
#endif

using std::cout;
using std::string;
using std::cin;


std::array<size_t, 4> get_nearest_four_ind(const std::vector<double>& dist_list, size_t dist_list_size, size_t centre);
double get_nearest_fifth_dist(const std::vector<double>& dist_list,size_t dist_list_size, size_t centre);
bool file_exists(const string& str);
double get_q_tet(const chemfiles::Frame& inputframe, const size_t centre, const std::array<size_t, 4>& nearest_four_indices, const std::vector<size_t>& oxy_ind_list);
double get_sk_val(const std::vector<double>& dist_list, size_t dist_list_size, size_t centre);
chemfiles::Vector3D get_centre_of_mass(const chemfiles::Frame& inputframe, const std::vector<double>& mass_list,const size_t natoms,const double mass_total);
bool is_within_bounds(const double value, const double low, const double high);
bool is_cell_at_edge(const std::vector<int>& neighbour_list);


constexpr double one_over_three = 1.0 / 3.0;
constexpr double three_over_eight = 3.0 / 8.0;
constexpr double num_to_dens = (18.01528 * 1.66054) ;

// Wrap everything into a try, except block to handle exceptions

int main(int argc, char** argv)
{
    // ready cin for throwing exceptions at bad input (as we are already using try block)
    cin.exceptions(std::ios_base::failbit);
try {
    // start counting time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // For task arg, only three options are possible -> OTO for Tetrahedral order, d5 for d5 parameter, Sk for translational tet. order
    const std::vector<string> task_allowed_args = { "OTO","d5","Sk","rhoV" };
    TCLAP::ValuesConstraint<string> allowed_args_constr(task_allowed_args);
    // Available command line arguments
    TCLAP::CmdLine cmdparser("Water order analysis", ' ', "1.0");
    TCLAP::ValueArg<size_t> end_frame_num_arg("", "stop", "Frame number to end calculation at; default is last frame (counting starts from 0)", false, 0, "positive integer", cmdparser);
    TCLAP::ValueArg<size_t> start_frame_num_arg("", "start", "Frame number to start calculation from; default is first frame (counting starts from 0)", false, 0, "positive integer", cmdparser);
    TCLAP::ValueArg<double> bin_size_arg("", "bin-width", "Size of each range considered for averaging", false, 0.2, "float (Angstrom)", cmdparser);
    TCLAP::ValueArg<double> rmax_arg("","rmax","Maximum distance from centre-of-mass considered for distribution (rmin = 0 always)",false,50.0,"float (Angstrom)",cmdparser);
    TCLAP::ValueArg<string> out_file_input_arg("o", "output-file", "Base name of output file; default is \'Water_order\'", false, "Water_order", "string", cmdparser);
    TCLAP::ValueArg<string> task_input_arg("t", "task", "Task requested to run: OTO = Oriental Tetrahedral Order parameter, d5 = d5 parameter, Sk = Translational Tetrahedral Order parameter, rhoV = Voronoi cell based water density; default is OTO", false, "OTO", &allowed_args_constr, cmdparser);
    TCLAP::ValueArg<string> o_name_input_arg("c", "oxygen-name", "Name of oxgyen atoms in the atom-info file", false, "NULL", "string", cmdparser);
    TCLAP::ValueArg<string> atinfo_input_arg("s", "atom-info", "Name of the file containing the atom information (topology)", false, "NULL", "string", cmdparser);
    TCLAP::ValueArg<string> trj_input_arg("f", "trajectory", "Name of the trajectory file", false, "NULL", "string", cmdparser);
    cmdparser.parse(argc, argv);


    /* Get the names of files containing atom information and the trajectory and the
        name of oxygen atoms. If not in command line, ask for them */
    string atinfo_file_name; // name of file with atom info
    string trj_file_name; // name of trajectory file
    string sel_string; // name of oxygen atom e.g. OH2
    string out_file_name = out_file_input_arg.getValue(); // name of the file for output
    double hist_end; // rmax value (i.e. end of histogram ranges)
    double bin_size; // bin-width in Angs
    

    // Check if no arguments are given

    if (atinfo_input_arg.getValue() == "NULL" && trj_input_arg.getValue() == "NULL") {
        cout << "No command line arguments for input files given, switching to keyboard input.\n";
        cout << "Please run " << argv[0] << " -h to see command line options for specifying input.\n";
        cout << "\n";
    }

    if (rmax_arg.isSet()) {
        hist_end = rmax_arg.getValue();
        cout << "Rmax read from command line: " << hist_end << " Angstrom\n";
    }
    else {
        cout << "Enter the maximum distance (Angstrom) from centre of mass that will be considered for outputting the data: ";
        cin >> hist_end;
    }

    if (bin_size_arg.isSet()) {
        bin_size = bin_size_arg.getValue();
        cout << "Width of each range read from command line: " << bin_size << " Angstrom\n";
    }
    else {
        cout << "Enter the size of each range for averaging: ";
        cin >> bin_size;
    }

    if (atinfo_input_arg.getValue() == "NULL") {
        cout << "Enter the name of file containing the atom information (i.e. .psf etc.): ";
        cin >> atinfo_file_name;
    }
    else {
        atinfo_file_name = atinfo_input_arg.getValue();
        cout << "Atom information (topology) file name read from command line: " << atinfo_file_name << "\n";
    }

    if (trj_input_arg.getValue() == "NULL") {
        cout << "Enter the name of the trajectory file (i.e. .dcd, .xtc etc.): ";
        cin >> trj_file_name;
    }
    else {
        trj_file_name = trj_input_arg.getValue();
        cout << "Trajectory file name read from command line: " << trj_file_name << "\n";
    }

    if (o_name_input_arg.getValue() == "NULL") {
        cout << "Enter the name of the oxygen atoms in the atom information file: ";
        cin >> sel_string;
    }
    else {
        sel_string = o_name_input_arg.getValue();
        cout << "Name of oxygen atom selected: " << sel_string << "\n";
    }

    if (out_file_input_arg.isSet()) {
        cout << "Output file name read from command line: " << out_file_name << "\n";
    }
    else {
        cout << "Using default output file name: " << out_file_name << "\n";
    }

    /* In this part, we read in the pdb or gro (or anything else) file that will
    give the names of atoms. This is crucial to detect the oxygen atoms */

    cout << "\n";
    cout << "Attempting to read the atom information file...\n"; // replace chemfiles with own code
    
    chemfiles::Trajectory atinfo(atinfo_file_name);
    chemfiles::Frame atinfo_frame = atinfo.read();
    if (atinfo_frame[0].name().empty()) {
        cout << "\n" << "Unable to read atom names from atom information file!" << "\n";
        cout << "File does not contain atom names, or is an unsupported format" << "\n" << "\n";
        cout << "Exiting" << "\n";
        exit(1);
    }
    else {
        cout << "  done" << "\n";
    }

    /*const std::vector<std::tuple<string, double>> ATOM_data = read_PSF(atinfo_file_name); // holds names and masses
    if (ATOM_data.size() == 0) {
        cout << "Unable to read atom names from atom information file!\n";
        exit(1);
    }
    else {
        cout << "...  done\n";
    } */

    // now get the atoms selected by sel_string and handle any errors
    cout << "Parsing selection string...";

    string sel_string_fin = "name " + sel_string;
    auto select_oxy = chemfiles::Selection(sel_string_fin);
    std::vector<size_t> oxy_matches = select_oxy.list(atinfo_frame); // vector contains indices of oxgyen atoms
    const auto oxy_size = oxy_matches.size(); // total number of oxygen atoms found

    /* std::vector <size_t> oxy_matches;
    for (size_t k = 0; k < ATOM_data.size(); ++k) {
        if (std::get<0>(ATOM_data.at(k))==sel_string) {
            oxy_matches.push_back(k);
        }
    }
    const auto oxy_size = oxy_matches.size(); */

    // if there are less than 4 O's then exit
    if (oxy_size < 5) {
        cout << "\n" << "There are 4 or less oxygen atoms (" << sel_string << ") in the system" << "\n";
        cout << "Order parameters can only be determined for >=5 oxygen atoms" << "\n";
        cout << "Exiting";
        exit(1);
    }
    else {
        cout << "  done" << "\n";
    }

    
     // Now read trajectory file
     cout << "Attempting to read the trajectory file..." << "\n";
     chemfiles::Trajectory trj(trj_file_name);
    if (trj.nsteps() <= 1) {
        cout << "\n" << "There is 1 or less frames in the trajectory file given!\n";
        cout << "Exiting\n";
        exit(1);
    }
    else {
        cout << "... done\n";
    }


    // check if the number of atoms is the same in atom info and trajectory files
    cout << "Checking compatibility of atom information and trajectory files...";
    // Testframe here!!
    trj.set_cell(chemfiles::UnitCell()); // ChemFiles has a problem with dcd files, set unit cell to no unit cell
    auto testframe = trj.read(); // for testing only!!
    if (testframe.size() == atinfo_frame.size()) {
        cout << "  done" << "\n";
    }
    else {
        cout << "\n" << "Number of atoms in the atom information file (" << atinfo_frame.size();
        cout << ") does not match the number of atoms in trajectory (" << testframe.size() << ") !\n";
        cout << "Please check the input files again.\n";
        cout << "Exiting" << std::endl;
        exit(1);
    }
    // now deal with the start and end frame
    size_t end_step;
    const size_t start_step = start_frame_num_arg.getValue();
    const size_t num_steps = trj.nsteps();

    if (start_frame_num_arg.isSet()) {
        if (start_step >= num_steps) {
            cout << "The value for first frame to be processed (" << start_step;
            cout << ") is larger than or equal to the total number of frames in the trajectory (" << num_steps << ").\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }
    }
    if (end_frame_num_arg.isSet()) {
        if (end_frame_num_arg.getValue() > num_steps) {
            cout << "The value for last frame to be processed (" << end_frame_num_arg.getValue();
            cout << ") is larger than the total number of frames in the trajectory (" << num_steps << ").\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }
        else {
            end_step = end_frame_num_arg.getValue();
        }
    }
    else {
        end_step = num_steps;
    }
    //const auto num_processed_frames = end_step - start_step;
    cout << "Calculation will start from frame " << start_step << " and will end at frame " << (end_step-1) << "\n";


    
    // make a list of masses
    // 
    auto natoms = atinfo_frame.size();
    std::vector<double> mass_list(natoms);
    for (size_t k = 0; k < natoms; ++k) {
        mass_list.at(k) = atinfo_frame[k].mass();
    }
    const double total_mass = std::accumulate(mass_list.begin(),mass_list.end(),0.0);
    // check if the file type actually gives mass. If it doesn't then zero mass would be set
    if (std::fabs(total_mass - 0.0) < 0.1) {
        cout << "The atom information / topology file does not contain masses of atoms!\n";
        cout << "Exiting";
        exit(1);
    }

    // set up histogramming

    const double hist_start = 0.0; // set rmin = 0 from COM
    if (hist_end <= 0.0) {
        cout << "Rmax value given (" << hist_end << ") is less than or equal to zero\n";
        cout << "Exiting";
        exit(1);
    }
    if (bin_size <= 0.0) {
        cout << "Error! bin_size given (" << bin_size << ") is less than or equal to zero\n";
        cout << "Exiting";
        exit(1);
    }
    size_t num_bins = (size_t)round((hist_end - hist_start) / bin_size);

    std::vector<double> bin_edges(num_bins + 1);
    // calculate bin edges
    for (int x = 0; x <= num_bins; ++x) {
        bin_edges.at(x) = hist_start + x * bin_size;
    }


    // This part runs if the Oriental Tetrahedral Order calculation is selected. The code for the 5th nearest
    // neighbour, d_5 is kept separate.
    if (task_input_arg.getValue() == "OTO") {
        cout << "Oriental Tetrahedral Order calculation requested \n";
        string out_file_name_oto = out_file_name + "_oto.csv";
        //check if the output file exists, if not then create it and open it
        if (file_exists(out_file_name_oto)) {
            // file exists, probably from previous run ... exit showing an error message
            cout << "Output file named " << out_file_name_oto << " exists!" << "\n";
            cout << "Past output files will not be overwritten." << "\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }
        std::ofstream outfile_oto;
        outfile_oto.open(out_file_name_oto.c_str());
        if (!outfile_oto) {
            cout << "File could not be opened!" << "\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }
        outfile_oto << std::fixed << std::showpoint; // fixed form and always output with decimal point
        outfile_oto << std::setprecision(6); // set 6 decimal points for float output
        // print the header into file
        outfile_oto << "r,qT\n";

        // for floating point overflow/underflow check
        std::feclearexcept(FE_UNDERFLOW);
        std::feclearexcept(FE_OVERFLOW);
        
        // Input checks done, now onto main calculations (OTO)

        std::vector<size_t> hist_counts(num_bins, 0); // histogram counts
        std::vector<double> hist_collect_q_tet(num_bins, 0.0); // collects the q_tet values        
        
        
        // iterate over whole trajectory
        #pragma omp parallel
        {
        if (omp_get_thread_num() == 0) {
            cout << "\nRunning on " << omp_get_num_threads() << " thread(s)\n\n";
        }
        std::vector<double> dist_oo(oxy_size); // vector for holding O-O distances
        std::array<size_t, 4> nearest_four_ind = { 0,0,0,0 }; //stores 4 nearest indices
        auto frame = chemfiles::Frame(); // holds the frame, necessary to have a dummy variable due to scoping
        double qval;
        std::vector<size_t> hist_counts_priv(num_bins, (size_t)0);
        std::vector<double> hist_collect_q_tet_priv(num_bins, 0.0);

        #pragma omp for schedule(static, 1)
        for (openmp_index_type n = start_step; n < end_step; ++n) {
            #pragma omp critical
            {
            frame = trj.read_step(n);
            }
            auto com_of_frame = get_centre_of_mass(frame, mass_list, natoms, total_mass); // gives COM for this frame
            auto positions_this_frame = frame.positions();
            // oxy_matches contains the actual indices of the oxygens, oxy_size contains the size
            for (size_t i = 0; i < oxy_size; ++i) {
                for (size_t j = 0; j < oxy_size; ++j) {
                    dist_oo[j] = frame.distance(oxy_matches[i], oxy_matches[j]); //gives distance between oxy_matches[i], oxy_matches[j] in Angs
                }
                nearest_four_ind = get_nearest_four_ind(dist_oo, oxy_size, i); //gives apparent indices back
                qval = get_q_tet(frame, i, nearest_four_ind, oxy_matches); // gives q_tet
                auto com_o_dist = (positions_this_frame[oxy_matches[i]] - com_of_frame).norm(); // get distance of O from COM
                for (int x = 0; x < num_bins; ++x) {
                    if (is_within_bounds(com_o_dist, bin_edges[x], bin_edges[x + 1])) {
                        hist_counts_priv[x]++;
                        hist_collect_q_tet_priv[x] += qval;
                        break;
                    }
                }
            }
            if (n % 200 == 0) {
                cout << "\rFrame " << n << " processed";
                //for (float count = 0.0f; count <= 1.0f; count += 0.01f) {
                //    if (n == (int)(num_processed_frames * count)) cout << (int)(count * 100) << "% done...";
                //}
            }
        }
        // sum up the histograms
        #pragma omp critical
        for (int x = 0; x < num_bins; ++x) {
            hist_counts[x] += hist_counts_priv[x];
            hist_collect_q_tet[x] += hist_collect_q_tet_priv[x];
        }
        }
        // parallel region processing ends here

        // check if there were overflow/underflow
        if (std::fetestexcept(FE_UNDERFLOW) || std::fetestexcept(FE_OVERFLOW)) {
            cout << "Floating point error!\nExiting";
            exit(1);
        }

        // write the histogram x and y into the file
        double y_intens = 0.0;
        for (int x = 0; x < num_bins; ++x) {
            if (hist_counts[x] == 0) {
                //y_intens = 0.0;
                outfile_oto << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << "\n";
            }
            else {
                y_intens = hist_collect_q_tet[x] / hist_counts[x];
                outfile_oto << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << y_intens << "\n";
            }
        }
        outfile_oto << std::endl; //flush the output stream
        outfile_oto.close();
        cout << "\nOriental Tetrahedral Order parameter calculation finished successfully!\n";
        cout << "Radial distribution written to " << out_file_name_oto << "\n";
    }
    
    // This part calculates the d5 parameter

    if (task_input_arg.getValue() == "d5") {
        cout << "d5 parameter calculation requested\n";
        string out_file_name_d5 = out_file_name + "_d5.csv";
        //check if the output file exists, if not then create it and open it
        if (file_exists(out_file_name_d5)) {
            // file exists, probably from previous run ... exit showing an error message
            cout << "Output file named " << out_file_name_d5 << " exists!" << "\n";
            cout << "Past output files will not be overwritten." << "\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }
        std::ofstream outfile_d5; // this is the output stream for the d5 output file
        outfile_d5.open(out_file_name_d5.c_str());
        if (!outfile_d5) {
            cout << "File could not be opened!" << "\n";
            cout << "Exiting" << std::endl;
        }
        outfile_d5 << std::fixed << std::showpoint; // fixed form and always output with decimal point
        outfile_d5 << std::setprecision(6); // set 6 decimal points for float output
        // print the header into file
        outfile_d5 << "r,d5\n";        

        // Input checks done, now onto main calculations (d5)

        // for floating point overflow/underflow check
        std::feclearexcept(FE_UNDERFLOW);
        std::feclearexcept(FE_OVERFLOW);
        
        // iterate over whole trajectory

        std::vector<size_t> hist_counts(num_bins, (size_t)0); // histogram counts
        std::vector<double> hist_collect_d5(num_bins,0.0); // collects the q_tet values
        
        #pragma omp parallel
        {
        std::vector<double> dist_oo(oxy_size); // holds distances between oxygens
        double d5_val; // holds d5 value
        auto frame = chemfiles::Frame(); // holds the frame, necessary to have a dummy variable due to scoping
        std::vector<size_t> hist_counts_priv(num_bins, (size_t)0);
        std::vector<double> hist_collect_d5_priv(num_bins, 0.0);
        if (omp_get_thread_num() == 0) {
            cout << "\nRunning on " << omp_get_num_threads() << " thread(s)\n\n";
        }

        #pragma omp for schedule(static, 1)
        for (openmp_index_type n = start_step; n < end_step; ++n) {
            #pragma omp critical
            {
            frame = trj.read_step(n);
            }
            auto com_of_frame = get_centre_of_mass(frame, mass_list, natoms, total_mass); // gives COM for this frame
            auto positions_this_frame = frame.positions();
            for (size_t i = 0; i < oxy_size; ++i) {
                for (size_t j = 0; j < oxy_size; ++j) {
                    dist_oo[j] = frame.distance(oxy_matches[i], oxy_matches[j]);
                }
                d5_val = get_nearest_fifth_dist(dist_oo, oxy_size, i); // holds d5 values
                auto com_o_dist = (positions_this_frame[oxy_matches[i]] - com_of_frame).norm(); // get distance of O from COM
                for (int x = 0; x < num_bins; ++x) {
                    if (is_within_bounds(com_o_dist, bin_edges[x], bin_edges[x + 1])) {
                        hist_counts_priv[x]++;
                        hist_collect_d5_priv[x] += d5_val;
                        break;
                    }
                }
            }
            if (n % 200 == 0) {
                cout << "\rFrame " << n << " processed";
                //for (float count = 0.0f; count <= 1.0f; count += 0.01f) {
                //    if (n == (short int)(num_processed_frames * count)) cout << (int)(count * 100) << "% done...";
                //}
            }
        }
        #pragma omp critical
        for (int x = 0; x < num_bins; ++x) {
            hist_counts[x] += hist_counts_priv[x];
            hist_collect_d5[x] += hist_collect_d5_priv[x];
        }
        }
        // parallel processing ends here

        // check if there were overflow/underflow
        if (std::fetestexcept(FE_UNDERFLOW) || std::fetestexcept(FE_OVERFLOW)) {
            cout << "Floating point error!\nExiting";
            exit(1);
        }

        // write the histogram x and y into the file
        double y_intens = 0.0;
        for (int x = 0; x < num_bins; ++x) {
            if (hist_counts[x] == 0) {
                //y_intens = 0.0;
                outfile_d5 << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << "\n";
            }
            else {
                y_intens = hist_collect_d5[x] / hist_counts[x];
                outfile_d5 << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << y_intens << "\n";
            }
        }
        outfile_d5 << std::endl; //flush the output stream
        outfile_d5.close();
        cout << "\nd5 parameter calculation finished successfully!\n";
        cout << "Radial distribution written to " << out_file_name_d5 << "\n";
    }

    if (task_input_arg.getValue() == "Sk") {
        cout << "Translational Tetrahedral Order parameter requested\n";
        string out_file_name_sk = out_file_name + "_Sk.csv";
        //check if the output file exists, if not then create it and open it
        if (file_exists(out_file_name_sk)) {
            // file exists, probably from previous run ... exit showing an error message
            cout << "Output file named " << out_file_name_sk << " exists!" << "\n";
            cout << "Past output files will not be overwritten." << "\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }
        std::ofstream outfile_sk; // this is the output stream for the d5 output file
        outfile_sk.open(out_file_name_sk.c_str());
        if (!outfile_sk) {
            cout << "File could not be opened!" << "\n";
            cout << "Exiting" << std::endl;
        }
        outfile_sk << std::fixed << std::showpoint; // fixed form and always output with decimal point
        outfile_sk << std::setprecision(6); // set 6 decimal points for float output
        // print the header into file
        outfile_sk << "r,Sk\n";
        // Input checks done, now onto main calculations (d5)

        // for floating point overflow/underflow check
        std::feclearexcept(FE_UNDERFLOW);
        std::feclearexcept(FE_OVERFLOW);

        // iterate over whole trajectory

        std::vector<size_t> hist_counts(num_bins,(size_t)0); // histogram counts
        std::vector<double> hist_collect_sk(num_bins,0.0); // collects the Sk values

        #pragma omp parallel
        {
        std::vector<double> dist_oo(oxy_size); //holds distances between oxygens
        double sk_val; //holds sk value for one oxygen
        auto frame = chemfiles::Frame();
        std::vector<size_t> hist_counts_priv(num_bins,(size_t)0);
        std::vector<double> hist_collect_sk_priv(num_bins,0.0);

        if (omp_get_thread_num() == 0) {
            cout << "\nRunning on " << omp_get_num_threads() << " thread(s)\n\n";
        }

        #pragma omp for schedule(static,1)
        for (openmp_index_type n = start_step; n < end_step; ++n) {
            #pragma omp critical
            {
            frame = trj.read_step(n);
            }
            auto com_of_frame = get_centre_of_mass(frame, mass_list, natoms, total_mass); // gives COM for this frame
            auto positions_this_frame = frame.positions();
            for (size_t i = 0; i < oxy_size; ++i) {
                for (size_t j = 0; j < oxy_size; ++j) {
                    dist_oo[j] = frame.distance(oxy_matches[i], oxy_matches[j]);
                }
                sk_val = get_sk_val(dist_oo,oxy_size,i);
                auto com_o_dist = (positions_this_frame[oxy_matches[i]] - com_of_frame).norm(); // get distance of O from COM
                for (int x = 0; x < num_bins; ++x) {
                    if (is_within_bounds(com_o_dist, bin_edges[x], bin_edges[x + 1])) {
                        hist_counts_priv[x]++;
                        hist_collect_sk_priv[x] += sk_val;
                        break;
                    }
                }
            }
            if (n % 200 == 0) {
                cout << "\rFrame " << n << " processed";
            }
        }
        // sum up histogram
        #pragma omp critical
        for (int x = 0; x < num_bins; ++x) {
            hist_counts[x] += hist_counts_priv[x];
            hist_collect_sk[x] += hist_collect_sk_priv[x];
        }
        }
        // parallel processing ends here
        // check if there were overflow/underflow
        if (std::fetestexcept(FE_UNDERFLOW) || std::fetestexcept(FE_OVERFLOW)) {
            cout << "Floating point error!\nExiting";
            exit(1);
        }
        // write the histogram x and y into the file
        double y_intens = 0.0;
        for (int x = 0; x < num_bins; ++x) {
            if (hist_counts[x] == 0) {
                //y_intens = 0.0;
                outfile_sk << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << "\n";
            }
            else {
                y_intens = hist_collect_sk[x] / hist_counts[x];
                outfile_sk << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << y_intens << "\n";
            }
        }
        outfile_sk << std::endl; //flush the output stream
        outfile_sk.close();
        cout << "\nSk parameter calculation finished successfully!\n";
        cout << "Radial distribution written to " << out_file_name_sk << "\n";

    }

    // Voronoi based density calculation
    if (task_input_arg.getValue() == "rhoV") {
        cout << "Voronoi cell based density requested\n";
        string out_file_name_rhov = out_file_name + "_rhoV.csv";
        //check if the output file exists, if not then create it and open it
        if (file_exists(out_file_name_rhov)) {
            // file exists, probably from previous run ... exit showing an error message
            cout << "Output file named " << out_file_name_rhov << " exists!" << "\n";
            cout << "Past output files will not be overwritten." << "\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }
        std::ofstream outfile_rhov; // this is the output stream for the d5 output file
        outfile_rhov.open(out_file_name_rhov.c_str());
        if (!outfile_rhov) {
            cout << "File could not be opened!" << "\n";
            cout << "Exiting" << std::endl;
        }
        outfile_rhov << std::fixed << std::showpoint; // fixed form and always output with decimal point
        outfile_rhov << std::setprecision(6); // set 6 decimal points for float output
        // print the header into file
        outfile_rhov << "r,rhoV\n";
        // Input checks done, now onto main calculations (Voronoi density)

        // for floating point overflow/underflow check
        std::feclearexcept(FE_UNDERFLOW);
        std::feclearexcept(FE_OVERFLOW);

        // we need all the heavy atoms, not just water for the Voronoi density
        std::vector<size_t> heavy_atom_matches;
        for (size_t k = 0; k < natoms; ++k) {
            if (mass_list.at(k) >= 3.0) { // heavy atom => anything heavier than hydrogen i.e. tritium?
                heavy_atom_matches.push_back(k);
            }
        }
        const auto n_heavy_atoms = heavy_atom_matches.size();
        if (n_heavy_atoms < oxy_size) {
            cout << "Number of heavy atoms is less that number of oxygen atoms.\n";
            cout << "Possibly an error in the masses in the topology.\n";
            cout << "Exiting" << std::endl;
            exit(1);
        }

        
        // iterate over whole trajectory
        std::vector<size_t> hist_counts(num_bins, 0); // histogram counts
        std::vector<double> hist_collect_rhov(num_bins); // collects the rho_V values

        // parallel part starts
        #pragma omp parallel
        {
        std::vector<double> rhov_val(n_heavy_atoms); // holds the rho_V value for all heavy atoms
        auto frame = chemfiles::Frame();
        double x, y, z; // required to hold the coordinates
        double x_max, y_max, z_max, x_min, y_min, z_min; // the max and min coordinates to get the box
        std::vector<size_t> hist_counts_priv(num_bins, 0); // thread private variable to hold data
        std::vector<double> hist_collect_rhov_priv(num_bins); // thread private variable

        if (omp_get_thread_num() == 0) {
            cout << "\nRunning on " << omp_get_num_threads() << " thread(s)\n\n";
        }
        #pragma omp for schedule(static,1)
        for (openmp_index_type n = start_step; n < end_step; ++n) {
            #pragma omp critical
            {
            frame = trj.read_step(n);
            }
            auto com_of_frame = get_centre_of_mass(frame, mass_list, natoms, total_mass); // gives COM for this frame
            auto positions_this_frame = frame.positions();
            x_max = positions_this_frame[0][0];
            y_max = positions_this_frame[0][1];
            z_max = positions_this_frame[0][2];
            x_min = x_max;
            y_min = y_max;
            z_min = z_max;

            // get the highest and lowest value of x,y and z for all heavy atoms
            for (size_t k = 0; k < n_heavy_atoms; ++k) {
                x = positions_this_frame[heavy_atom_matches[k]][0];
                y = positions_this_frame[heavy_atom_matches[k]][1];
                z = positions_this_frame[heavy_atom_matches[k]][2];
                if (x > x_max) x_max = x;
                if (y > y_max) y_max = y;
                if (z > z_max) z_max = z;
                if (x < x_min) x_min = x;
                if (y < y_min) y_min = y;
                if (z < z_min) z_min = z;
            }
            const auto span_x = x_max - x_min; // these hold the span in each axes
            const auto span_y = y_max - y_min;
            const auto span_z = z_max - z_min;
            unsigned int block_size = (unsigned int)std::cbrt(5.0 * span_x * span_y * span_z / (double)n_heavy_atoms); // ideally each block should have 5 particles
            if (block_size==0) {
                block_size = 1;
            }
            // a box big enough to store all the particles
            voro::container mybox(x_min - 0.5*span_x, x_max + 0.5*span_x, y_min - 0.5*span_y, y_max + 0.5*span_y, z_min - 0.5*span_z, z_max + 0.5*span_z, block_size, block_size, block_size, false, false, false, 1);
            for (size_t k = 0; k < n_heavy_atoms; ++k) {
                x = positions_this_frame[heavy_atom_matches[k]][0];
                y = positions_this_frame[heavy_atom_matches[k]][1];
                z = positions_this_frame[heavy_atom_matches[k]][2];
                mybox.put(k, x, y, z);
            }
            voro::c_loop_all loop_cursor(mybox);
            voro::voronoicell_neighbor this_cell;
            std::vector<int> neighbour_list;
            if (loop_cursor.start()) do if (mybox.compute_cell(this_cell, loop_cursor)) {
                auto k = loop_cursor.pid();
                this_cell.neighbors(neighbour_list);
                if (is_cell_at_edge(neighbour_list)) {
                    rhov_val.at(k) = 0.0;
                }
                else {
                    rhov_val.at(k) = num_to_dens * 1.0 / this_cell.volume(); // volume is in number per Angstrom^3
                }
            } while (loop_cursor.inc());
            // now deal with not all heavy atoms being oxygens
            for (size_t k = 0; k < n_heavy_atoms; ++k) {
                // find if the index in heavy atom match is also in oxygen match => it is an oxygen
                if (std::find(oxy_matches.begin(), oxy_matches.end(), heavy_atom_matches[k]) != oxy_matches.end()) {
                    auto com_o_dist = (positions_this_frame[heavy_atom_matches[k]] - com_of_frame).norm(); // get distance of O from COM
                    // populate the histogram
                    for (int q = 0; q < num_bins; ++q) {
                        if (is_within_bounds(com_o_dist, bin_edges[q], bin_edges[q + 1])) {
                            hist_counts_priv[q]++;
                            hist_collect_rhov_priv[q] += rhov_val[k];
                            break;
                        }
                    }
                    
                }
            }
            if (n % 200 == 0) {
                cout << "\rFrame " << n << " processed";
            }
        }
        // sum up private histograms
        #pragma omp critical
        for (int q = 0; q < num_bins; ++q) { // x is reserved for coordinates
            hist_counts[q] += hist_counts_priv[q];
            hist_collect_rhov[q] += hist_collect_rhov_priv[q];
        }
        }
        // parallel part ends
        if (std::fetestexcept(FE_UNDERFLOW) || std::fetestexcept(FE_OVERFLOW)) {
            cout << "Floating point error!\nExiting";
            exit(1);
        }
        // write the histogram x and y into the file
        double y_intens = 0.0;
        for (int x = 0; x < num_bins; ++x) {
            if (hist_counts[x] == 0) {
                //y_intens = 0.0;
                outfile_rhov << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << "\n";
            }
            else {
                y_intens = hist_collect_rhov[x] / hist_counts[x];
                outfile_rhov << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << y_intens << "\n";
            }
        }
        outfile_rhov << std::endl; //flush the output stream
        outfile_rhov.close();
        cout << "\nVoronoi cell based density parameter calculation finished successfully!\n";
        cout << "Radial distribution written to " << out_file_name_rhov << "\n";


    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    cout << std::setprecision(4);
    cout << "Total wall time elapsed (sec) = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())/1000000.0 << std::endl;
    // execution finished!!

} catch (const chemfiles::Error &chemex) {
    cout << "\n" << "Error in chemfiles:" << chemex.what() << std::endl;
    exit(1);
} catch (TCLAP::ArgException &tclex) {
    cout << "\n" << "Error in TCLAP:" << tclex.error() << "for arg " << tclex.argId() << std::endl;
    exit(1);
} catch (std::ios_base::failure &ioex) {
    cout << "\n" << "Error in keyboard input:" << ioex.what() << "\n";
    cout << "Please check that input data is of correct type!"<< std::endl;
    exit(1);
} catch (...) {
    cout << "\n" << "Other exception caught" << std::endl;
    exit(1);
}
    return 0;
}


/// <summary>
/// Get the lowest five distances from the dist_list, because the lowest would of
/// course, be zero. 0th number of the array would always be zero.The function
/// actually returns the apparent indices of the nearest four oxygen, *not* the distances themselves.
/// </summary>
/// <param name="dist_list">- Vector list of distances between oxygens, using their apparent indices</param>
/// <param name="dist_list_size">- Size of dist_list</param>
/// <param name="centre">- apparent index of central oxygen</param>
/// <returns>(array of 4 integers) the 4 nearest integers </returns>
std::array<size_t, 4> get_nearest_four_ind(const std::vector<double>& dist_list, size_t dist_list_size, size_t centre) {

    // the code is sort of hacky, but since only 4 nearest are required, it works
    // a grid based search would be faster but is much more difficult to write

    std::array<size_t, 4> near_ind = { 0,0,0,0 }; // these hold the apparent indices of atoms with minimum distances
    // indices should be like [ind0], ind1, ind2, ind3, ind4
    // initialize the min values to -1 to check if they have been assigned or not
    double min1=-1.0, min2=-1.0, min3=-1.0, min4=-1.0; // these hold the actual distances

    // we know min0 == 0, min0 is not required
    // first get ind0 i.e. centre -> obtained from function parameter

    // get first minimum value after 0
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre) {
            if (min1 == -1.0 || dist_list[m] < min1) {
                min1 = dist_list[m];
                near_ind[0] = m;
            }
        }
    }
    // get second minimum value
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != near_ind[0]) {
            if (min2 == -1.0 || dist_list[m] < min2) {
                min2 = dist_list[m];
                near_ind[1] = m;
            }
        }
    }
    // get third minimum value
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != near_ind[0] && m != near_ind[1]) {
            if (min3 == -1.0 || dist_list[m] < min3) {
                min3 = dist_list[m];
                near_ind[2] = m;
            }
        }
    }
    // get fourth minimum value
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != near_ind[0] && m != near_ind[1] && m != near_ind[2]) {
            if (min4 == -1.0 || dist_list[m] < min4) {
                min4 = dist_list[m];
                near_ind[3] = m;
            }
        }
    }



    // simple checking -> not necessary to check all combinations
    // just a defense against wrong code
    if (centre == near_ind[0] || near_ind[0] == near_ind[1] || centre == near_ind[0] || near_ind[0] == near_ind[2] || near_ind[1] == near_ind[2]) {
        cout << "Error in finding four nearest neighbours" << "\n";
        cout << "Exiting...";
        exit(1);
    }

    return near_ind; // return the indices
}

/// <summary>
/// This function returns the distance to the 5th nearest neighbour of the selected (central) oxygen
/// </summary>
/// <param name="dist_list">Vector list of distances between oxygens, using their apparent indices</param>
/// <param name="dist_list_size">Size of dist_list</param>
/// <param name="centre">apparent index of central oxygen</param>
/// <returns>(double) the value of d5 parameter for central oxygen</returns>
double get_nearest_fifth_dist(const std::vector<double>& dist_list, size_t dist_list_size, size_t centre)
{
    size_t ind1 = 0, ind2 = 0, ind3 = 0, ind4 = 0, ind5 = 0;
    double min1 = -1.0, min2 = -1.0, min3 = -1.0, min4 = -1.0, min5 = -1.0;

    // we know min0 == 0, and ind0 is centre

    // get ind1, min1
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre) {
            if (min1 == -1.0 || dist_list[m] < min1) {
                min1 = dist_list[m];
                ind1 = m;
            }
        }
    }
    // get ind2, min2
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != ind1) {
            if (min2 == -1.0 || dist_list[m] < min2) {
                min2 = dist_list[m];
                ind2 = m;
            }
        }
    }
    // get ind3, min3
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != ind1 && m != ind2) {
            if (min3 == -1.0 || dist_list[m] < min3) {
                min3 = dist_list[m];
                ind3 = m;
            }
        }
    }
    // get ind4, min4
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != ind1 && m != ind2 && m != ind3) {
            if (min4 == -1.0 || dist_list[m] < min4) {
                min4 = dist_list[m];
                ind4 = m;
            }
        }
    }
    // get ind5, min5
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != ind1 && m != ind2 && m != ind3 && m != ind4) {
            if (min5 == -1.0 || dist_list[m] < min5) {
                min5 = dist_list[m];
                ind5 = m;
            }
        }
    }

    return min5;
}



/// <summary>
/// This function checks if a file exists or not, by attempting to read the file with 
/// std::ifstream. If the file exists, return true.
/// </summary>
/// <param name="str">String containing the file name</param>
/// <returns>(bool) true if file exists</returns>
bool file_exists(const string& str)
{
    std::ifstream fs(str.c_str());
    return fs.is_open();
}

/// <summary>
/// This function evaluates the tetrahedral order parameter.
/// </summary>
/// <param name="inputframe">Passes a reference to the current frame (because we need the frame.angle() method)</param>
/// <param name="centre">Apparent index of the central atom</param>
/// <param name="nearest_four_indices">Array of the apparent indices of four nearest</param>
/// <param name="oxy_ind_list">Vector that converts apparent indices to actual indices</param>
/// <returns>(double) value of qT for central oxygen</returns>
double get_q_tet(const chemfiles::Frame& inputframe, const size_t centre, const std::array<size_t, 4>& nearest_four_indices, const std::vector<size_t>& oxy_ind_list) {
    //          O0
    //          |
    //    O2----O(centre)----O1
    //          |
    //          O3
    // 
    double q_tet_collect = 0.0;
    for (short int j = 0; j < 3; ++j) {
        for (short int k = j + 1; k < 4; ++k) {
            q_tet_collect += std::pow(
            std::cos(
            inputframe.angle(oxy_ind_list[nearest_four_indices[j]], oxy_ind_list[centre], oxy_ind_list[nearest_four_indices[k]])
            )
            + one_over_three, 2
            );
            //cout << j << " " << centre << k << "\n";
        }
    }
    //cout << q_tet_collect;

    return (1.0 - (three_over_eight * q_tet_collect));
}
/// <summary>
/// A function that determines the Sk (translational tetrahedral order parameter)
/// </summary>
/// <param name="dist_list">- a vector containing the distances from central oxygen to all other oxygens</param>
/// <param name="dist_list_size">- size of dist_list</param>
/// <param name="centre">- apparent index of central atom</param>
/// <returns>- (double) the Sk value for central oxygen</returns>
double get_sk_val(const std::vector<double>& dist_list, size_t dist_list_size, size_t centre)
{
    double sk_collect = 0.0;
    std::array<size_t, 4> near_ind = { 0,0,0,0 }; // these hold the apparent indices of atoms with minimum distances
    // indices should be like [ind0], ind1, ind2, ind3, ind4
    // initialize the min values to -1 to check if they have been assigned or not
    double min1=-1.0, min2=-1.0, min3=-1.0, min4=-1.0; // these hold the actual distances

    // we know min0 == 0, min0 is not required
    // first get ind0 i.e. centre -> obtained from function parameter

    // get first minimum value after 0
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre) {
            if (min1 == -1.0 || dist_list[m] < min1) {
                min1 = dist_list[m];
                near_ind[0] = m;
            }
        }
    }
    // get second minimum value
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != near_ind[0]) {
            if (min2 == -1.0 || dist_list[m] < min2) {
                min2 = dist_list[m];
                near_ind[1] = m;
            }
        }
    }
    // get third minimum value
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != near_ind[0] && m != near_ind[1]) {
            if (min3 == -1.0 || dist_list[m] < min3) {
                min3 = dist_list[m];
                near_ind[2] = m;
            }
        }
    }
    // get fourth minimum value
    for (size_t m = 0; m < dist_list_size; ++m) {
        if (m != centre && m != near_ind[0] && m != near_ind[1] && m != near_ind[2]) {
            if (min4 == -1.0 || dist_list[m] < min4) {
                min4 = dist_list[m];
                near_ind[3] = m;
            }
        }
    }

    // simple checking -> not necessary to check all combinations
    // just a defense against wrong code
    if (centre == near_ind[0] || near_ind[0] == near_ind[1] || centre == near_ind[0] || near_ind[0] == near_ind[2] || near_ind[1] == near_ind[2]) {
        cout << "Error in finding four nearest neighbours" << "\n";
        cout << "Exiting...";
        exit(1);
    }
    auto r_avg = (min1 + min2 + min3 + min4) / 4;
    std::array<double, 4> r_k = { {min1,min2,min3,min4} };
    for (short int z = 0; z < 4; ++z) {
        sk_collect += std::pow((r_k[z] - r_avg), 2); // calculate sum
    }
    sk_collect = sk_collect / (4 * std::pow(r_avg, 2));

    return (1.0 - (one_over_three * sk_collect)); // return the Sk
}

/* This function evaluates the centre of mass of a frame obtained from ChemFiles
    \param inputframe - reference to the ChemFiles frame object
    \param mass_list - a vector containing the masses of all atoms, from topology
    \param natoms - total number of atoms in the frame
    \param mass_total - total mass of the atoms in the frame */
chemfiles::Vector3D get_centre_of_mass(const chemfiles::Frame& inputframe, const std::vector<double>& mass_list, const size_t natoms, const double mass_total)
{
    auto frame_com = chemfiles::Vector3D();
    auto atom_positions = inputframe.positions();
    for (size_t i = 0; i < natoms; ++i) {
        frame_com += atom_positions[i] * mass_list[i];
    }
    frame_com = frame_com / mass_total;
    return frame_com;
}

/* This function evaluates if a number is within a given range [low,high)
    \param value - the specified number
    \param low - lower end of range
    \param high - higher end of range */
bool is_within_bounds(const double value, const double low, const double high)
{
    return !(value < low) && (value < high);
}
/* This function evaluates if the cell is at the edge of the container box. If it at the edge, one of the neighbouring 
 particles would have ID of -1 to -6 (to indicate the 6 walls)
    \param neighbour_list - A std::vector<int> containing information about the neighbouring particles */
bool is_cell_at_edge(const std::vector<int>& neighbour_list)
{
    if (neighbour_list.size() == 0) {
        cout << "Neighbour list has 0 size, exiting.\n"<<std::endl;
        exit(1);
    }
    for (const auto z : neighbour_list) {
        if (z < 0) {
            return true;
        }
        else {
            return false;
        }
    }
}

