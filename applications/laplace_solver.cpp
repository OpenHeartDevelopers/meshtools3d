#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include "Chrono.hpp"
#include "GetPot.hpp"
#include "LaplacePipeline.hpp"
#include "Mesh.hpp"

// Collect all values that immediately follow a given flag in argv.
// Handles repeated flags: --zero-bc a.vtx --zero-bc b.vtx
static std::vector<std::string> collectFlagValues(int argc, char** argv,
                                                   const std::string& flag)
{
    std::vector<std::string> results;
    for (int i = 1; i < argc - 1; ++i) {
        if (std::string(argv[i]) == flag)
            results.push_back(std::string(argv[i + 1]));
    }
    return results;
}

int main(int argc, char** argv)
{
    GetPot command_line(argc, argv);

    if (command_line.search(2, "-h", "--help") || argc == 1) {
        std::cout << "laplace_solver\n"
                  << "Runs the Laplace solver (and optionally thickness evaluation)\n"
                  << "on an existing CARP mesh.\n\n"
                  << "Usage:\n"
                  << "  ./laplace_solver --zero-bc <vtx> --one-bc <vtx>\n"
                  << "                   -mesh_dir <dir> -mesh_name <name>\n"
                  << "                   -out_dir <dir>  -out_name <name>\n"
                  << "                   [-f <data_file>] [options]\n\n"
                  << "BC flags (can be repeated for multiple node sets):\n"
                  << "  --zero-bc <vtx_file>           node set assigned value 0.0\n"
                  << "  --one-bc  <vtx_file>           node set assigned value 1.0\n\n"
                  << "Mesh / output:\n"
                  << "  -mesh_dir <dir>                directory containing CARP mesh\n"
                  << "  -mesh_name <name>              CARP mesh basename\n"
                  << "  -out_dir <dir>                 output directory\n"
                  << "  -out_name <name>               output file basename\n\n"
                  << "Options:\n"
                  << "  -f <data_file>                 GetPot file for [laplacesolver] params\n"
                  << "  --no-thickness                 solve only, skip thickness evaluation\n"
                  << "  --swap-regions                 swap endo/epi propagation direction\n"
                  << "  --thickness-algorithm N, -thickalgo N\n"
                  << "                                 1=Bishop (default), 2=Corrado\n"
                  << "  --vtk                          write VTK output\n"
                  << "  --vtk-binary                   binary VTK (requires --vtk)\n"
                  << "  --potential                    write potential output\n"
                  << "  -v, --verbose                  verbose solver output\n";
        return 0;
    }

    // Optional data file for solver params only
    std::string data_file_name = command_line.follow("", 2, "-f", "--file");

    std::string mesh_dir  = command_line.follow(".", 2, "-mesh_dir",  "--mesh_directory");
    std::string mesh_name = command_line.follow("mesh", 2, "-mesh_name", "--mesh_basename");
    std::string out_dir   = command_line.follow(".", 2, "-out_dir",  "--output_directory");
    std::string out_name  = command_line.follow("laplace_out", 2, "-out_name", "--output_name");

    // BCs
    LaplaceBCConfig bc;
    bc.zero_bc_files = collectFlagValues(argc, argv, "--zero-bc");
    bc.one_bc_files  = collectFlagValues(argc, argv, "--one-bc");
    bc.swap_regions  = command_line.search(1, "--swap-regions");

    if (bc.zero_bc_files.empty() || bc.one_bc_files.empty()) {
        std::cerr << "ERROR: at least one --zero-bc and one --one-bc file are required.\n";
        return 1;
    }

    // Solver params — from data file if provided, otherwise defaults
    LaplaceParams params;
    if (!data_file_name.empty()) {
        GetPot param_file(data_file_name.c_str());
        params.abs_toll    = param_file("laplacesolver/abs_toll",    1e-6);
        params.rel_toll    = param_file("laplacesolver/rel_toll",    1e-6);
        params.itr_max     = param_file("laplacesolver/itr_max",     500);
        params.dimKrilovSp = param_file("laplacesolver/dimKrilovSp", 150);
        params.verbose     = param_file("laplacesolver/verbose",     0);
    }
    if (command_line.search(2, "-v", "--verbose"))
        params.verbose = 1;

    // Output config
    OutputConfig out_cfg;
    out_cfg.out_dir        = out_dir;
    out_cfg.out_name       = out_name;
    out_cfg.eval_thickness = !command_line.search(1, "--no-thickness");
    out_cfg.out_vtk        = command_line.search(1, "--vtk");
    out_cfg.out_vtk_binary = command_line.search(1, "--vtk-binary");
    out_cfg.out_potential  = command_line.search(1, "--potential");

    int algo = 1;
    if (command_line.search(2, "-thickalgo", "--thickness-algorithm"))
        algo = command_line.follow(1, 2, "-thickalgo", "--thickness-algorithm");
    out_cfg.thickness_algorithm = algo;

    // Create output directory
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (system(mkdir_cmd.c_str()) != 0) {
        std::cerr << "ERROR: could not create output directory: " << out_dir << "\n";
        return 1;
    }

    // Load mesh
    std::string mesh_path = mesh_dir + "/" + mesh_name;
    std::cout << "Reading mesh: " << mesh_path << "\n";
    Chrono chrono;
    chrono.start();
    Mesh mesh(mesh_path);
    chrono.stop();
    std::cout << "Mesh read in " << chrono << "\n";

    runLaplacePipeline(mesh, bc, params, out_cfg);

    return 0;
}
