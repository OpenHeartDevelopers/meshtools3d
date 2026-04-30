// parfile_builder — emit a meshtools3d parameter file with sensible defaults.
//
// Pure stdlib. No CGAL, no m3d/, no third-party dependencies. The schema
// (sections, keys, defaults) tracks docs/parameter_file_schema.md (currently
// .claude/m3d_python_params.md) verbatim so the C++ writer and the future
// Python writer in pycemrg-meshing stay in lockstep.

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace {

using KV      = std::pair<std::string, std::string>;
using Section = std::pair<std::string, std::vector<KV>>;
using Schema  = std::vector<Section>;

const Schema& defaults()
{
    static const Schema s = {
        {"segmentation", {
            {"seg_dir",                "./"},
            {"seg_name",               "seg_final_smooth_corrected.inr"},
            {"mesh_from_segmentation", "1"},
            {"boundary_relabeling",    "0"},
        }},
        {"meshing", {
            {"facet_angle",         "30"},
            {"facet_size",          "0.8"},
            {"facet_distance",      "4"},
            {"cell_rad_edge_ratio", "2.0"},
            {"cell_size",           "0.8"},
            {"rescaleFactor",       "1000"},
        }},
        {"laplacesolver", {
            {"abs_toll",    "1e-6"},
            {"rel_toll",    "1e-6"},
            {"itr_max",     "700"},
            {"dimKrilovSp", "500"},
            {"verbose",     "1"},
        }},
        {"others", {
            {"eval_thickness", "0"},
        }},
        {"output", {
            {"outdir",          "./myocardium_OUT"},
            {"name",            "heart_mesh"},
            {"out_medit",       "0"},
            {"out_carp",        "1"},
            {"out_carp_binary", "0"},
            {"out_vtk",         "1"},
            {"out_vtk_binary",  "0"},
            {"out_potential",   "0"},
        }},
    };
    return s;
}

void printHelp(std::ostream& os)
{
    os <<
"Usage: parfile_builder [OPTIONS]\n"
"\n"
"Write a meshtools3d parameter file with sensible defaults.\n"
"Output goes to stdout unless -o/--out is given.\n"
"\n"
"Options:\n"
"  -o, --out PATH            Output file path\n"
"  --seg_dir DIR             Set segmentation.seg_dir\n"
"  --seg_name FILE           Set segmentation.seg_name\n"
"  --out_dir DIR             Set output.outdir\n"
"  --out_name NAME           Set output.name\n"
"  --set SECTION.KEY=VAL     Override any key (repeatable)\n"
"  -h, --help                Show this help\n"
"\n"
"Examples:\n"
"  parfile_builder > heart.par\n"
"  parfile_builder -o heart.par --seg_dir /data --out_dir /data/out\n"
"  parfile_builder --set meshing.facet_size=0.5 --set output.out_vtk=1\n";
}

bool applyOverride(Schema& schema,
                   const std::string& section,
                   const std::string& key,
                   const std::string& value)
{
    for (auto& [sname, kvs] : schema)
    {
        if (sname != section) continue;
        for (auto& kv : kvs)
        {
            if (kv.first == key)
            {
                kv.second = value;
                return true;
            }
        }
        std::cerr << "error: unknown key '" << key
                  << "' in section [" << section << "]\n"
                  << "  valid keys:";
        for (const auto& kv : kvs) std::cerr << " " << kv.first;
        std::cerr << "\n";
        return false;
    }
    std::cerr << "error: unknown section '" << section << "'\n"
              << "  valid sections:";
    for (const auto& [sname, _] : schema) std::cerr << " " << sname;
    std::cerr << "\n";
    return false;
}

bool parseSetArg(const std::string& arg,
                 std::string& section,
                 std::string& key,
                 std::string& value)
{
    auto dot = arg.find('.');
    auto eq  = arg.find('=');
    if (dot == std::string::npos || eq == std::string::npos || dot >= eq)
    {
        std::cerr << "error: --set expects SECTION.KEY=VALUE, got '" << arg << "'\n";
        return false;
    }
    section = arg.substr(0, dot);
    key     = arg.substr(dot + 1, eq - dot - 1);
    value   = arg.substr(eq + 1);
    return true;
}

void writeFile(std::ostream& os, const Schema& schema)
{
    bool first = true;
    for (const auto& [section, kvs] : schema)
    {
        if (!first) os << "\n";
        first = false;
        os << "[" << section << "]\n";
        for (const auto& [k, v] : kvs)
            os << k << " = " << v << "\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    Schema schema = defaults();
    std::string outPath;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        auto needArg = [&](const char* flag) -> std::string {
            if (i + 1 >= argc)
            {
                std::cerr << "error: " << flag << " requires a value\n";
                std::exit(1);
            }
            return argv[++i];
        };

        if (arg == "-h" || arg == "--help")
        {
            printHelp(std::cout);
            return 0;
        }
        else if (arg == "-o" || arg == "--out")
        {
            outPath = needArg("-o/--out");
        }
        else if (arg == "--seg_dir")
        {
            if (!applyOverride(schema, "segmentation", "seg_dir",  needArg("--seg_dir")))  return 1;
        }
        else if (arg == "--seg_name")
        {
            if (!applyOverride(schema, "segmentation", "seg_name", needArg("--seg_name"))) return 1;
        }
        else if (arg == "--out_dir")
        {
            if (!applyOverride(schema, "output", "outdir",         needArg("--out_dir")))  return 1;
        }
        else if (arg == "--out_name")
        {
            if (!applyOverride(schema, "output", "name",           needArg("--out_name"))) return 1;
        }
        else if (arg == "--set")
        {
            std::string s, k, v;
            if (!parseSetArg(needArg("--set"), s, k, v)) return 1;
            if (!applyOverride(schema, s, k, v))         return 1;
        }
        else
        {
            std::cerr << "error: unknown argument '" << arg << "'\n"
                      << "Run with --help for usage.\n";
            return 1;
        }
    }

    if (outPath.empty())
    {
        writeFile(std::cout, schema);
    }
    else
    {
        std::ofstream f(outPath);
        if (!f)
        {
            std::cerr << "error: could not open '" << outPath << "' for writing\n";
            return 1;
        }
        writeFile(f, schema);
    }
    return 0;
}
