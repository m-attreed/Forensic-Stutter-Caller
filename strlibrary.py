import math
import csv
import os, sys
import unittest


loci = ["AMEL", "D3S1358", "D1S1656", "D2S441", "D10S1248",
        "D13S317", "Penta E", "D16S539", "D18S51", "D2S1338", "CSF1PO",
        "Penta D", "TH01", "vWA", "D21S11", "D7S820", "D5S818",
        "TPOX", "DYS391", "D8S1179", "D12S391", "D19S433", "FGA", "D22S1045"]

channel_dictionary = {
        "Sample Name": "", "AMEL": "Blue", "D3S1358": "Blue",
        "D1S1656": "Blue", "D2S441": "Blue", "D10S1248": "Blue",
        "D13S317": "Blue", "Penta E": "Blue", "D16S539": "Green",
        "D18S51": "Green", "D2S1338": "Green", "CSF1PO": "Green",
        "Penta D": "Green", "TH01": "Yellow", "vWA": "Yellow",
        "D21S11": "Yellow", "D7S820": "Yellow", "D5S818": "Yellow",
        "TPOX": "Yellow", "DYS391": "Yellow", "D8S1179": "Red",
        "D12S391": "Red", "D19S433": "Red", "FGA": "Red", "D22S1045": "Red"}


class AlleleUnit:

    def __init__(self, allele, locus=None):
        if allele in ["X", "Y", "INC", "OB", "OL"]:
            self.locusType = "Character"
            self.allele = allele
        else:
            self.locusType = "Number"
            bps, reps = math.modf(float(allele))
            self.basepairs = round(bps * 10)
            self.repeats = int(reps)
            self.basepairsPerRepeat = self.assign_bp_per_repeat(locus)

        self.locus = locus

    def assign_bp_per_repeat(self, locus):
        if locus is not None:
            return self.look_up_bp_in_repeat(locus)
        else:
            return 0

    nonTetramerDict = {"Penta D": 5, "Penta E": 5, "D22S1045": 3}

    def confirm_locus(self, locus):
        return (locus in loci)
    
    def look_up_bp_in_repeat(self, locus):
        if locus in self.nonTetramerDict:
            return self.nonTetramerDict[locus]
        else:
            return 4

    def convert_to_bp(self):
        return (self.repeats * self.basepairsPerRepeat) + self.basepairs

    def convert_to_repeats(self, repeats):
        return float(str(repeats//self.basepairsPerRepeat)+"."+str(repeats%self.basepairsPerRepeat))

    def __eq__(self, other):
        return self.__str__() == other.__str__()

    def __ne__(self, other):
        return (not self.__eq__(other))

    def __lt__(self, other):
        if self.locusType == "Character" or other.locusType == "Character":
            return self.__str__() < other.__str__()
        else:
            return float(self.__repr__()) < float(other.__repr__())

    def __hash__(self):
        return hash(self.__repr__())

    def __repr__(self):
        if self.locusType == "Number":
            if self.basepairs != 0:
                allele = str(self.repeats)+"."+str(self.basepairs)
            else:
                allele = str(self.repeats)
        else:
            allele = self.allele
        return allele

    def __str__(self):
        if self.locusType == "Number":
            if self.basepairs != 0:
                allele = "Allele: "+str(self.repeats)+"."+str(self.basepairs)
            else:
                allele = "Allele: "+str(self.repeats)
        else:
            allele = "Allele: "+self.allele
        return allele

    def add(self, other):
        basepairs1 = self.convert_to_bp()
        basepairs2 = other.convert_to_bp()
        totalBasepairs = basepairs1 + basepairs2
        self.repeats = totalBasepairs//self.basepairsPerRepeat
        self.basepairs = totalBasepairs%self.basepairsPerRepeat

    def subtract(self, other):
        basepairs1 = self.convert_to_bp()
        basepairs2 = other.convert_to_bp()
        totalBasepairs = basepairs1 - basepairs2
        self.repeats = totalBasepairs//self.basepairsPerRepeat
        self.basepairs = totalBasepairs%self.basepairsPerRepeat



class Profile:
    """
    The Profile class expects a list of strings.
    The first string is the sample name
    and the rest are the loci as strings.
    It expects comma seperated alleles.
    It accepts profiles as Genotypes or Phenotypes
    it internally converts them to phenotypes.
    """

    def __init__(self, profilesArray, mixName=None):
        """
        Work on this section I am trying to make it so that I can take multiple input profiles
        or just one profile and if multiple combine them, then process the profiles as usual
        if multiple it should get a mixName and apply that to the mixture

        mixName should be the same as the mix name in the sample name So that it can be looked up
        as in does mixName appear in sample name in case multiple mixtures are found in a single
        plate of data
        """

        # I want to throw an error if there are multiple profiles but no mixName
        self.mixName = mixName
        tempProfileArray = []
        if len(profilesArray) > 1:
            tempProfileArray = list(map(",".join,zip(*profilesArray)))
        else:
            tempProfileArray = profilesArray[0]

        self.sampleName, profileArray = tempProfileArray[0], tempProfileArray[1:]
        if mixName != None:
            self.sampleName = self.mixName

        profileArray = [profile.split(',') for profile in profileArray]

        profileArray = list(zip(loci, profileArray))
        alleleArray = [[AlleleUnit(allele, locus[0]) for allele in locus[1]]for locus in profileArray]
        alleleArrayToPhenotype = []

        for locus in alleleArray:
            tempLocus = list(set(locus))
            tempLocus.sort()
            alleleArrayToPhenotype.append(tempLocus)

        self.profile = dict(zip(loci, alleleArrayToPhenotype))

    def combineWith(self, other):
        for value in self.profile:
            tempValuesSelf = self.profile[value]
            tempValuesOther = other.profile[value]
            tempValuesSelf.extend(tempValuesOther)

            tempValuesSelf = list(set(tempValuesSelf))
            tempValuesSelf.sort()
            if AlleleUnit("INC") in tempValuesSelf and len(tempValuesSelf) > 1:
                tempValuesSelf.remove(AlleleUnit("INC"))
            self.profile[value] = tempValuesSelf
        return self.profile

    def __str__(self):
        header = "Sample Name\t"+"\t".join(loci)
        profilePrint = self.sampleName+"\t"
        for locus in loci:
            alleles = self.profile[locus]
            if len(alleles) > 1:
                alleles = [allele.__repr__() for allele in alleles]
                profilePrint = profilePrint+",".join(alleles)
                profilePrint = profilePrint+"\t"
            else:
                profilePrint = profilePrint+alleles[0].__repr__()+"\t"
        return header+"\n"+profilePrint+"\n"


class ProfileDB:
    def __init__(self, profileCSVFile):
        tempProfilesList = []
        with open(profileCSVFile, newline='') as data:
            data_reader = csv.reader(data, delimiter='\t')
            for line in data_reader:
                tempProfilesList.append(line)
        self.profilesHeaderRow = tempProfilesList.pop(0)
        tempProfilesDB = []
        for profile in tempProfilesList:
            tempProfile = Profile([profile])
            tempProfilesDB.append(tempProfile)
        self.profilesDB = tempProfilesDB

    def getProfile(self, searchName):
        foundProfile = None
        for profile in self.profilesDB:
            if profile != None and profile.sampleName == searchName:
                foundProfile = profile
        return foundProfile

    def addMixes(self, mixFile):
        """This takes a list of profiles and turns in into a single profile
        that combines all the profiles; profiles are 25 elements"""

        tempMixArray = []

        with open(mixFile, newline='') as data:
            data_reader = csv.reader(data, delimiter='\t')
            for line in data_reader:
                tempMixArray.append(line)
        numProfiles = len(tempMixArray)



        for x in range(1, numProfiles):
            tempMixProfile = None
            tempMixName = tempMixArray[x][0]+"-"+tempMixArray[x][1]
            for y in range(2, 2 + int(tempMixArray[x][1])):
                if tempMixProfile == None:
                    tempMixProfile = self.getProfile(tempMixArray[x][y])
                    tempMixProfile.sampleName = tempMixName
                else:
                    combinedProfile = tempMixProfile.combineWith(self.getProfile(tempMixArray[x][y]))
                    tempMixProfile.profile = combinedProfile

            if tempMixProfile != None:
                self.profilesDB.append(tempMixProfile)

    def __str__(self):
        headerRow = "\t".join(self.profilesHeaderRow)+"\n"
        profileDBString = headerRow
        for profile in self.profilesDB:
            tempProfileString = profile.__str__().split("\n")[1]
            profileDBString += tempProfileString+"\n"
        return profileDBString


class ReportDB:
    """ follow the ProfileDB class and pull from the stutter flagger file
    this should make up dictionary of file names that have a list made up
    of all the lines that have that file name
    """

    def __init__(self, file):
        inputFile = []
        self.fileName = file

        with open(self.fileName, newline='') as data:
            data_reader = csv.reader(data, delimiter='\t')
            for line in data_reader:
                inputFile.append(line)

        tempHeaderRow = inputFile.pop(0)
        tempHeaderRow[0] = "#"
        tempHeaderRow.append("NOC")
        tempHeaderRow.append("Type")


        self.reportHeaderRow = tempHeaderRow


        self.report = inputFile

        self.Marker = self.reportHeaderRow.index("Marker")
        self.Dye = self.reportHeaderRow.index("Dye")
        self.Size = self.reportHeaderRow.index("Size")
        self.Allele = self.reportHeaderRow.index("Allele")
        self.Sample_Comments = self.reportHeaderRow.index("Sample Comments")
        self.Height = self.reportHeaderRow.index("Height")
        self.NOC = self.reportHeaderRow.index("NOC")
        self.Program_Output = self.reportHeaderRow.index("Type")
        self.profilesInDB = self.profile_codes()
        self.sampleProperties = self.define_sample_properties()

        self.sampleList = self.sample_set()
        self.samplesSorted = self.collect_sample_data()
        self.samplePropertiesDict = self.make_properties_dict()

        self.samplePullupDict = self.make_pullup_dict()

    def profile_codes(self):
        """Takes a list of profiles as lines from the input and splits them
        using the underscore character and takes the sample name. This excludes
        the ladder and the Amp Neg samples. All other samples including the
        Amp Pos sample are added to the set."""

        """Changed to use only manually specified profiles for mixtures 
        where samples are not in the file name anymore."""
        profilesSet = set()
        for sample in self.report:
            profileCode = sample[1].split("_")[1]
            if profileCode not in ['Allelic Ladder', 'Amp Neg']:
                profilesSet.add(profileCode)
        return profilesSet

    def sample_set(self):
        """This produces a list of all the samples in the input file
        by file name."""
        sampleSet = set()
        for line in self.report:
            sampleSet.add(line[1])
        return sampleSet

    def collect_sample_data(self):
        """This divides the input file into lists by sample file name.
        This is used to feed the data sample by sample through the application."""
        dataBySample = []

        for sample in self.sampleList:
            sampleDataOnly = []
            for line in self.report:
                if line[1] == sample:
                    sampleDataOnly.append(line)
            dataBySample.append(sampleDataOnly)
        return dataBySample

    def make_properties_dict(self):
        propertiesDict = {}

        for sample in self.sampleList:
            propertiesDict[sample] = SampleProperties(sample)
        return propertiesDict

    def make_pullup_dict(self):
        pullupDict = {}

        for sample in self.sampleList:
            pullupDict[sample] = SamplePullup(sample)
        return pullupDict


    def mark_parent_peaks(self, profilesDB):
        """This function marks the parent peaks for the current data set.
        This iterates through the lines of the input and marks parent peaks.

        It skips the Amelogenin and DYS391 loci.

        It checks for the following issues: overlapping, dropout,
        saturation, and ILS failure. While determining if there is dropout it
        sets the flags if dropout is found."""
        for sampleSet in self.samplesSorted:
            sampleName = sampleSet[1][1]
            NOC = 1
            if "Mix" in sampleSet[1][1].split("_")[1]:
                name = sampleSet[1][1].split("_")[1]
                if len(sampleSet[1][1].split("_")[2].split("-")) == 1:
                    NOC = len(sampleSet[1][1].split("_")[3].split("-"))
                else:
                    NOC = len(sampleSet[1][1].split("_")[2].split("-"))
                sampleForData = name+"-"+str(NOC)
            elif "Amp_Pos" in sampleSet[1][1]:
                sampleForData = "Amp Pos"
            elif "Positive" in sampleSet[1][1]:
                sampleForData = "Amp Pos"
            elif "Amp_Neg" in sampleSet[1][1]:
                sampleForData = "Amp Neg"
            elif "Ladder" in sampleSet[1][1]:
                sampleForData = "Ladder"
            else:
                sampleForData = sampleSet[1][1].split("_")[1]
            profileForData = profilesDB.getProfile(sampleForData)
            for peak in sampleSet:

                peak.append(str(NOC))
                if peak[self.Marker] != '' and sampleForData not in ["Ladder", "Amp Neg"] \
                    and peak[self.Sample_Comments] not in \
                        ["ILS Failure", "ILS Fails", "Misplating Fails", "Size Call Failed"]:
                    currentAllele = AlleleUnit(peak[self.Allele])

                    if currentAllele in profileForData.profile[peak[self.Marker]]:
                        peak.append("Par")
                        if peak[self.Dye] == "Blue":
                            self.samplePullupDict[sampleName].Blue.append(peak[self.Size])
                        elif peak[self.Dye] == "Green":
                            self.samplePullupDict[sampleName].Green.append(peak[self.Size])
                        elif peak[self.Dye] == "Yellow":
                            self.samplePullupDict[sampleName].Yellow.append(peak[self.Size])
                        elif peak[self.Dye] == "Red":
                            self.samplePullupDict[sampleName].Red.append(peak[self.Size])

                        currentPropDict = self.samplePropertiesDict[sampleName].loci[peak[self.Marker]]
                        #self.samplePropertiesDict[sampleName].loci[peak[self.Marker]].Peak_BP.append(peak[self.Size])
                        currentPropDict.Peak_BP.append(peak[self.Size])
                        #self.samplePropertiesDict[sampleName].loci[peak[self.Marker]].Peak_Profiles.append(peak[self.Allele])
                        currentPropDict.Peak_Profiles.append(peak[self.Allele])
                        #self.samplePropertiesDict[sampleName].loci[peak[self.Marker]].Channel = peak[self.Channel]
                        currentPropDict.Channel = peak[self.Dye]

                    else:
                        peak.append("X")

                elif "Fail" in peak[self.Sample_Comments]:
                    peak.append("Fail")
                else:
                    peak.append("X")


    # Rework the following three functions to work with the current structure of the program
    def in_stutter_position(self, parent, position, allele):
        """This function determines if a non-parent peak is in stutter position
        based on allele bin."""
        locus = allele[self.Marker]
        if allele[self.Allele] not in ["OL"]:
            return AlleleUnit(parent,locus).add(AlleleUnit(position,locus)) == AlleleUnit(allele[self.Allele])
        else:
            return False

    def is_loci_of_interest(self, locus):
        """Determines if the locus is one of the loci we are analyzing in this
        study only for tetramer locations only. Penta D and E and D22S1045 are
        handled separately."""
        return locus not in ['', "AMEL"]

    def is_allele_markable(self, flagData, peak, peakNumber):
        """Determines if the current potential stutter peak is markable
        based on the following criteria."""
        return peakNumber in flagData[peak[self.Marker]]["Called Peaks"] \
               and peak[self.Program_Output] == "X"

    def is_within_bp_range(self, parentPeak, peak, stutterPos):
        """Determines if the potential stutter peak is in stutter position
        based on the size of the parent peaks at the locus."""
        return ((float(parentPeak) - 0.5) + stutterPos) \
               <= float(peak[self.Size]) \
               <= ((float(parentPeak) + 0.5) + stutterPos)


    def mark_stutter(self, profilesDB):
        """
        This function marks all stutter peaks.
        Peaks are labeled db b hb or f followed by the number of the parent peak.

        remember that some locations aren't repeats of 4bp
        Penta loci have 5 basepair repeats and D22S1045 has 3 basepair repeats.
        These loci are separated from the tetramer loci.

        This function checks flags to determine if the allele can be used.
        """

        for sampleSet in self.samplesSorted:
            #NOC = -1
            if "Mix" in sampleSet[1][1].split("_")[1]:
                name = sampleSet[1][1].split("_")[1]
                NOC = len(sampleSet[1][1].split("_")[2].split("-"))
                sampleForData = name + "-" + str(NOC)
            elif "Amp_Pos" in sampleSet[1][1]:
                sampleForData = "Amp Pos"
                #NOC = 1
            elif "Ladder" in sampleSet[1][1]:
                sampleForData = "Ladder"
            else:
                sampleForData = sampleSet[1][1].split("_")[1]
                #NOC = 1
            profileForData = profilesDB.getProfile(sampleForData)

            for peak in sampleSet:
                #peak[self.NOC] = NOC
                if self.is_loci_of_interest(peak[self.Marker]):
                    alleles = self.samplePropertiesDict[sampleSet[1][1]].loci[peak[self.Marker]].Peak_Profiles
                    peakSizes = self.samplePropertiesDict[sampleSet[1][1]].loci[peak[self.Marker]].Peak_BP
                    for x in range(len(alleles)):
                        repeatMultiple = 0
                        if peak[self.Marker] not in ["D22S1045", "Penta D", "Penta E"]:
                            repeatMultiple = 4
                        elif peak[self.Marker] == "D22S1045":
                            repeatMultiple = 3
                        else:
                            repeatMultiple = 5
                        parentPeak = peakSizes[x]
                        parentAllele = alleles[x]



                        if self.is_within_bp_range(parentPeak, peak, (repeatMultiple * -1)) \
                                or self.in_stutter_position(parentAllele, -1, peak):
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "b"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output]+",b"

                        elif self.is_within_bp_range(parentPeak, peak, (repeatMultiple * -2)) \
                                or self.in_stutter_position(parentAllele, -2, peak):
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "db"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output] + ",db"

                        elif repeatMultiple == 4 and (self.is_within_bp_range(parentPeak, peak, (repeatMultiple * -0.5)) \
                                or self.in_stutter_position(parentAllele, -0.2, peak)):
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "hb"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output] + ",hb"

                        elif self.is_within_bp_range(parentPeak, peak, (repeatMultiple * 1)) \
                                or self.in_stutter_position(parentAllele, 1, peak):
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "f"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output] + ",f"

    def mark_pullup(self):
        """
        This function takes the pullupList from definePullup and use a test
        """
        for sampleSet in self.samplesSorted:

            sampleName = sampleSet[1][1]


            for peak in sampleSet:
                if peak[self.Dye] == "Blue":
                    for peakSize in self.samplePullupDict[sampleName].Green + self.samplePullupDict[sampleName].Yellow + self.samplePullupDict[sampleName].Red:
                        if float(peak[self.Size]) - 0.5 <= float(peakSize) \
                                <= float(peak[self.Size]) + 0.5:
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "pullup"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output] + ",pullup"
                elif peak[self.Dye] == "Green":
                    for peakSize in self.samplePullupDict[sampleName].Blue + self.samplePullupDict[sampleName].Yellow + self.samplePullupDict[sampleName].Red:
                        if float(peak[self.Size]) - 0.5 <= float(peakSize) \
                                <= float(peak[self.Size]) + 0.5:
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "pullup"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output] + ",pullup"
                elif peak[self.Dye] == "Yellow":
                    for peakSize in self.samplePullupDict[sampleName].Blue + self.samplePullupDict[sampleName].Green + self.samplePullupDict[sampleName].Red:
                        if float(peak[self.Size]) - 0.5 <= float(peakSize) \
                                <= float(peak[self.Size]) + 0.5:
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "pullup"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output] + ",pullup"
                elif peak[self.Dye] == "Red":
                    for peakSize in self.samplePullupDict[sampleName].Blue + self.samplePullupDict[sampleName].Green + self.samplePullupDict[sampleName].Yellow:
                        if float(peak[self.Size]) - 0.5 <= float(peakSize) \
                                <= float(peak[self.Size]) + 0.5:
                            if peak[self.Program_Output] == "X":
                                peak[self.Program_Output] = "pullup"
                            else:
                                peak[self.Program_Output] = peak[self.Program_Output] + ",pullup"


    def write_output(self):
        outputFileName = self.fileName.rsplit('.', 1)[0]

        outputFile = open(outputFileName + "_newoutput.tsv", "w")
        outputFile.write("\t".join(self.reportHeaderRow)+"\n")
        for line in self.report:
            outputFile.write("\t".join(line)+"\n")
        outputFile.close()


    def define_sample_properties(self):
        return {0: 0}


    def __str__(self):
        headerRow = "\t".join(self.reportHeaderRow)+"\n"
        reportDBString = headerRow
        for line in self.report:
            tempReportString = "\t".join(line)
            reportDBString += tempReportString+"\n"
        return reportDBString


class LocusProperties:
    def __init__(self):

        self.Channel = ""
        self.Saturation = False
        self.Drop_Out = False
        self.Peak_BP = []
        self.Peak_Profiles = []

    def __str__(self):
        return "Channel: "+self.Channel+"\n"+"Saturation: "+str(self.Saturation)+"\n"+\
               "Drop out: "+str(self.Drop_Out)+"\n"+"Peak Profiles: "+str(self.Peak_Profiles)+"\n"+\
               "Peak BP: "+str(self.Peak_BP)+"\n"


class SampleProperties:
    def __init__(self, sampleName):
        self.sampleName = sampleName

        self.loci = {
            "AMEL": LocusProperties(), "D3S1358": LocusProperties(), "D1S1656": LocusProperties(),
            "D2S441": LocusProperties(), "D10S1248": LocusProperties(), "D13S317": LocusProperties(),
            "Penta E": LocusProperties(), "D16S539": LocusProperties(), "D18S51": LocusProperties(),
            "D2S1338": LocusProperties(), "CSF1PO": LocusProperties(), "Penta D": LocusProperties(),
            "TH01": LocusProperties(), "vWA": LocusProperties(), "D21S11": LocusProperties(),
            "D7S820": LocusProperties(), "D5S818": LocusProperties(), "TPOX": LocusProperties(),
            "DYS391": LocusProperties(), "D8S1179": LocusProperties(), "D12S391": LocusProperties(),
            "D19S433": LocusProperties(), "FGA": LocusProperties(), "D22S1045": LocusProperties()}

    def __str__(self):


        outputString = ""
        for locus in self.loci:
            outputString += locus
            outputString += "\n"
            outputString += self.loci[locus].__str__()


        return "Sample name: "+self.sampleName+"\n"+outputString

class SamplePullup:
    """
    This function takes the flagdata at the "Peak BP" location for
    all keys in dictionary and makes a set of those data points in
    case any overlap this list is passed to the next function to filter
    the all called peaks
    """

    def __init__(self, sampleName):
        self.sampleName = sampleName
        self.Blue = []
        self.Green = []
        self.Yellow = []
        self.Red = []

    def __str__(self):
        blue = "\, ".join(self.Blue)+"\n"
        green = "\, ".join(self.Green)+"\n"
        yellow = "\, ".join(self.Yellow)+"\n"
        red = "\, ".join(self.Red)+"\n"
        return "Pullup BP: "+"\n"+blue+green+yellow+red


class TestAlleleUnit(unittest.TestCase):

    def test_microvariant(self):
        self.assertFalse(AlleleUnit(14.1, "TPOX") == AlleleUnit(14, "TPOX"))
        self.assertFalse(AlleleUnit(11.3, "TPOX") == AlleleUnit(11.2, "TPOX"))
        self.assertTrue(AlleleUnit(10, "TPOX") == AlleleUnit(10.0, "TPOX"))


if __name__ == "__main__":
    unittest.main()
