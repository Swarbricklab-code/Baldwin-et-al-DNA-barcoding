##ClonTracer count barcodes

#Script and algorithm sourced from Addgene 
#https://www.addgene.org/pooled-library/clontracer/
#originally cited in Bhang et al Nat Med 2015

# Script for the analysis of DNA barcode sequencing. 
# Input files: fastq
# Output: .csv files
# Script run with all default parameters using python 2



__author__ = 'krishvi7'
import sys
import getopt
from modules import classes

clte_object = classes.ClonTracerCountBarcodesExperiment()


def main(argv):
    arg_single = ":".join(clte_object.args_dict.values())+":" + "".join(clte_object.args_dict_paramless.values())
    arg_verbose = map(lambda x: x + "=", clte_object.args_dict.keys())
    arg_verbose.extend(map(lambda x: x, clte_object.args_dict_paramless.keys()))

    try:
        opts, args = getopt.getopt(argv, arg_single, arg_verbose)
    except getopt.GetoptError as exp:
        classes.Common.print_error("Invalid input argument: "+exp.msg)

    if len(opts) == 0:
        clte_object.usage()

    for opt, arg in opts:
        if opt in map(lambda x: "-"+x, clte_object.args_dict.values()) or \
           opt in map(lambda x: "--"+x, clte_object.args_dict.keys()):
            if opt in map(lambda x: "-"+x, clte_object.args_dict.values()):
                clte_member = clte_object.args_dict.keys()[clte_object.args_dict.values().index(opt[1:])]
            else:
                clte_member = opt[2:]

            try:
                tryout = classes.Common().dataType_dict[clte_member](arg)
            except ValueError as val_err:
                classes.Common.print_error("Invalid argument "+arg+" for parameter "+opt)

            if clte_member in classes.Common().range_dict.keys():
                if classes.Common().dataType_dict[clte_member](arg) not in classes.Common().range_dict[clte_member]:
                    raise classes.Common.print_error("Argument for "+opt+" must be " +
                                                         ",".join(map(str, classes.Common().range_dict[clte_member])))
            setattr(clte_object, clte_member, classes.Common().dataType_dict[clte_member](arg))

        elif opt in map(lambda x: "-"+x, clte_object.args_dict_paramless.values()) or \
           opt in map(lambda x: "--"+x, clte_object.args_dict_paramless.keys()):
            if opt in map(lambda x: "-"+x, clte_object.args_dict_paramless.values()):
                clte_member = clte_object.args_dict_paramless.keys()[clte_object.args_dict_paramless.values().index(opt[1:])]
            else:
                clte_member = opt[2:]

            setattr(clte_object, clte_member, classes.Common().paramless_dict[clte_member])
        else:
            classes.Common.print_error("Invalid input argument "+opt)

    try:
        print("Evaluating arguments")
        clte_object.check_input_directory()
        clte_object.check_output_directory()
        clte_object.check_known_contaminant_file()
        clte_object.count_barcodes_on_samples()

        print("Suppress quality summary: "+str(clte_object.suppress_quality_cw_summary))
        print("Library mode: "+str(clte_object.library_mode))

        print("")
        if not clte_object.suppress_quality_cw_summary:
            clte_object.cw_quality_summary()
        if not (clte_object.suppress_tall_skinny or clte_object.library_mode):
            clte_object.tall_skinny()
        if not (clte_object.suppress_matrix or clte_object.library_mode):
            clte_object.matricize()

    except:
        raise

if __name__ == "__main__":
    main(sys.argv[1:])


