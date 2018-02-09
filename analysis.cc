#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <getopt.h>
#include <string>

#include <hloop.h>

#include "core.h"

int main(int argc, char **argv)
{
    TROOT Analysis("Analysis","compiled analysis macro");
    gStyle->SetOptStat(0);

    // argc is the number of arguments in char* array argv
    // CAUTION: argv[0] contains the progname
    // argc has to be nargs+1

    /* Flag set by ‘--verbose’. */
    AnaParameters anapars;
    anapars.verbose_flag = 0;
    anapars.nomdc_flag = 0;
    anapars.nosys_flag = 0;
    anapars.norpc_flag = 0;
    anapars.decay_only_flag = 0;

    anapars.events = -1;

    anapars.config = "";
    anapars.outfile = "output.root";

    int c;
    while (1)
    {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            {"verbose",   no_argument,       &(anapars.verbose_flag), 1},
            {"brief",     no_argument,       &(anapars.verbose_flag), 0},
            {"nomdc",     no_argument,       &(anapars.nomdc_flag), 1},
            {"nosys",     no_argument,       &(anapars.nosys_flag), 1},
            {"norpc",     no_argument,       &(anapars.norpc_flag), 1},
            {"decay",     no_argument,       &(anapars.decay_only_flag), 1},
            /* These options don’t set a flag.
             *              We distinguish them by their indices. */
            {"config",    required_argument, 0, 'c'},
            {"events",    required_argument, 0, 'e'},
            {"output",    required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "c:e:o:",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 'c':
                anapars.config = optarg;
                break;

            case 'e':
                anapars.events = atol(optarg);
                break;

            case 'o':
                anapars.outfile = optarg;
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort ();
        }
    }

    /* Instead of reporting ‘--verbose’
     *      and ‘--brief’ as they are encountered,
     *           we report the final status resulting from them. */
    if (anapars.verbose_flag)
        puts ("verbose flag is set");

    HLoop * loop = new HLoop(kTRUE);

    /* Print any remaining command line arguments (not options). */
    if (optind < argc)
    {
        Bool_t ret = kFALSE;
        while (optind < argc)
        {
            TString infile = argv[optind++];
            if        (infile.Contains(","))        ret = loop->addMultFiles(infile);
            else if (infile.Contains(".root"))    ret = loop->addFiles(infile);
            else                                ret = loop->addFilesList(infile);

            if (!ret)
            {
                std::cerr << "READBACK: ERROR : cannot find inputfiles : " << infile.Data() << endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

    core(loop, anapars);

    exit(0);
}


