// Water_order.cpp : This file contains the 'main' function. Program execution begins and ends there.
/* Oriental tetrahedral order parameter calculator
 Written by Shoubhik Raj Maiti, Aug 2021 */

/* Reading of chemical structure files is done with Chemfiles by Guillaume Fraux
 Command line argument parsing is done with TCLAP by Mike Smoot and Daniel Aarno */

// 

#include <chemfiles.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <tclap/CmdLine.h>
#include <cmath>
#include <fstream>
#include <utility>
#include <tuple>
#include <iomanip>
#include <numeric>
#include <omp.h>
#include <cfenv>

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

std::vector<std::tuple<string, double>> read_PSF(string infile_name);
std::array<size_t, 4> get_nearest_four_ind(const std::vector<double>& dist_list, size_t dist_list_size, size_t centre);
double get_nearest_fifth_dist(const std::vector<double>& dist_list,size_t dist_list_size, size_t centre);
bool file_exists(const string& str);
double get_q_tet(const chemfiles::Frame& inputframe, const size_t centre, const std::array<size_t, 4>& nearest_four_indices, const std::vector<size_t>& oxy_ind_list);
chemfiles::Vector3D get_centre_of_mass(const chemfiles::Frame& inputframe, const std::vector<double>& mass_list,const size_t natoms,const double mass_total);
bool is_within_bounds(const double value, const double low, const double high);


constexpr double one_over_three = 1.0 / 3.0;
constexpr double three_over_eight = 3.0 / 8.0;

// Wrap everything into a try, except block to handle exceptions

int main(int argc, char** argv)
{
    // ready cin for throwing exceptions at bad input (as we are already using try block)
    cin.exceptions(std::ios_base::failbit);
try {
    // For task arg, only two options are possible -> OTO for Tetrahedral order and d5 for d5 parameter
    const std::vector<string> task_allowed_args = { "OTO","d5" };
    TCLAP::ValuesConstraint<string> allowed_args_constr(task_allowed_args);
    // Available command line arguments
    TCLAP::CmdLine cmdparser("Water order analysis", ' ', "1.0");
    TCLAP::ValueArg<size_t> end_frame_num_arg("", "stop", "Frame number to end calculation at; default is last frame (counting starts from 0)", false, 0, "positive integer", cmdparser);
    TCLAP::ValueArg<size_t> start_frame_num_arg("", "start", "Frame number to start calculation from; default is first frame (counting starts from 0)", false, 0, "positive integer", cmdparser);
    TCLAP::ValueArg<double> bin_size_arg("", "bin-width", "Size of each range considered for averaging", false, 0.2, "float (Angstrom)", cmdparser);
    TCLAP::ValueArg<double> rmax_arg("","rmax","Maximum distance from centre-of-mass considered for distribution (rmin = 0 always)",false,50.0,"float (Angstrom)",cmdparser);
    TCLAP::ValueArg<string> out_file_input_arg("o", "output-file", "Base name of output file; default is \'Water_order\'", false, "Water_order", "string", cmdparser);
    TCLAP::ValueArg<string> task_input_arg("t", "task", "Task requested to run: OTO = Oriental Tetrahedral Order, d5 = d5 parameter; default is OTO", false, "OTO", &allowed_args_constr, cmdparser);
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
    
    /* chemfiles::Trajectory atinfo(atinfo_file_name);
    chemfiles::Frame atinfo_frame = atinfo.read();
    if (atinfo_frame[0].name().empty()) {
        cout << "\n" << "Unable to read atom names from atom information file!" << "\n";
        cout << "File does not contain atom names, or is an unsupported format" << "\n" << "\n";
        cout << "Exiting" << "\n";
        exit(1);
    }
    else {
        cout << "  done" << "\n";
    } */

    const std::vector<std::tuple<string, double>> ATOM_data = read_PSF(atinfo_file_name); // holds names and masses
    if (ATOM_data.size() == 0) {
        cout << "Unable to read atom names from atom information file!\n";
        exit(1);
    }
    else {
        cout << "...  done\n";
    }

    // now get the atoms selected by sel_string and handle any errors
    cout << "Parsing selection string...";

    /* string sel_string_fin = "name " + sel_string;
    auto select_oxy = chemfiles::Selection(sel_string_fin);
    std::vector<size_t> oxy_matches = select_oxy.list(atinfo_frame); // vector contains indices of oxgyen atoms
    const auto oxy_size = oxy_matches.size(); // total number of oxygen atoms found */
    std::vector <size_t> oxy_matches;
    for (size_t k = 0; k < ATOM_data.size(); ++k) {
        if (std::get<0>(ATOM_data.at(k))=="OH2") {
            oxy_matches.push_back(k);
        }
    }
    const auto oxy_size = oxy_matches.size();

    // if there are less than 4 O's then exit
    if (oxy_size < 5) {
        cout << "\n" << "There are 4 or less oxygen atoms (" << sel_string << ") in the system" << "\n";
        cout << "Tetrahedral order parameter can only be determined for >=5 oxygen atoms" << "\n";
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
    if (testframe.size() == ATOM_data.size()) {
        cout << "  done" << "\n";
    }
    else {
        cout << "\n" << "Number of atoms in the atom information file (" << ATOM_data.size();
        cout << ") does not match the number of atoms in trajectory (" << trj.read().size() << ") !\n";
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
    const auto num_processed_frames = end_step - start_step;
    cout << "Calculation will start from frame " << start_step << " and will end at frame " << end_step-1 << "\n";


    
    // make a list of masses
    // 
    auto natoms = ATOM_data.size();
    std::vector<double> mass_list(natoms);
    for (size_t k = 0; k < natoms; ++k) {
        mass_list.at(k) = std::get<1>(ATOM_data.at(k));
    }
    const double total_mass = std::accumulate(mass_list.begin(),mass_list.end(),0.0);

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
    int num_bins = (int)round((hist_end - hist_start) / bin_size);

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
        std::vector<double> hist_collect_q_tet(num_bins); // collects the q_tet values        
        
        
        // iterate over whole trajectory
        #pragma omp parallel
        {
        if (omp_get_thread_num() == 0) {
            cout << "\nRunning on " << omp_get_num_threads() << " threads\n\n";
        }
        std::vector<double> dist_oo(oxy_size); // vector for holding O-O distances
        std::array<size_t, 4> nearest_four_ind = { 0,0,0,0 }; //stores 4 nearest indices
        auto frame = chemfiles::Frame(); // holds the frame, necessary to have a dummy variable due to scoping
        double qval;

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
                        #pragma omp atomic
                        hist_counts[x]++;
                        #pragma omp atomic
                        hist_collect_q_tet[x] += qval;
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
                y_intens = 0.0;
            }
            else {
                y_intens = hist_collect_q_tet[x] / hist_counts[x];
            }
            outfile_oto << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << y_intens << "\n";
        }
        outfile_oto << std::endl; //flush the output stream
        outfile_oto.close();
        cout << "\nOTO calculation finished successfully!\n";
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

        std::vector<size_t> hist_counts(num_bins, 0); // histogram counts
        std::vector<double> hist_collect_d5(num_bins); // collects the q_tet values
        
        #pragma omp parallel
        {
        std::vector<double> dist_oo(oxy_size); // holds distances between oxygens
        double d5_val; // holds d5 value
        auto frame = chemfiles::Frame(); // holds the frame, necessary to have a dummy variable due to scoping
        if (omp_get_thread_num() == 0) {
            cout << "\nRunning on " << omp_get_num_threads() << " threads\n\n";
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
                        #pragma omp atomic
                        hist_counts[x]++;
                        #pragma omp atomic
                        hist_collect_d5[x] += d5_val;
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
                y_intens = 0.0;
            }
            else {
                y_intens = hist_collect_d5[x] / hist_counts[x];
            }
            outfile_d5 << (bin_edges[x] + bin_edges[x + 1]) / 2 << "," << y_intens << "\n";
        }
        outfile_d5 << std::endl; //flush the output stream
        outfile_d5.close();
        cout << "\nd5 parameter calculation finished successfully!\n";
        cout << "Radial distribution written to " << out_file_name_d5 << "\n";
    }
    
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

/* Get the lowest five distances from the dist_list, because the lowest would of
 course, be zero. 0th number of the array would always be zero. The function
 actually returns the apparent indices of the nearest four oxygen, *not* the distances themselves.
    \param dist_list - Vector list of distances between oxygens, using their apparent indices
    \param dist_list_size - Size of dist_list
    \param centre - apparent index of central oxygen  */
std::array<size_t, 4> get_nearest_four_ind(const std::vector<double>& dist_list, size_t dist_list_size, size_t centre) {

    // the code is sort of hacky, but since only 4 nearest are required, it works

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

/* This function returns the distance to the 5th nearest neighbour of the selected (central) oxygen
    \param dist_list - Vector list of distances between oxygens, using their apparent indices
    \param dist_list_size - Size of dist_list
    \param centre - apparent index of central oxygen */
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


/* This function checks if a file exists or not, by attempting to read the file with
std::ifstream. If the file exists, return true.
    \param str - String containing the file name */
bool file_exists(const string& str)
{
    std::ifstream fs(str.c_str());
    return fs.is_open();
}

/* This function evaluates the tetrahedral order parameter.
    \param inputframe - Passes a reference to the current frame (because we need the frame.angle() method)
    \param centre - Apparent index of the central atom
    \param nearest_four_indices - Array of the apparent indices of four nearest
    \param oxy_ind_list - Vector that converts apparent indices to actual indices */
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


bool is_within_bounds(const double value, const double low, const double high)
{
    return !(value < low) && (value < high);
}

