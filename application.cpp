/*---------------------------------------------------------------------
@brief Project 1: DNA Profiling
@author Uvaish Bakaliya
@date Sep 15 2022

@brief Program Components:
    load_db
    load_dna
    display
    process
    search
    load_multi_dnas (own component)


Creative component:

    Name of the component: 'load_multi_dnas'

    Instructions:
        Entering the beginning and ending numbers in the range of those reads
        those DNAs files and pushes the data of those files into a DNA data
        container that can be used later to display and can also be used in the
        many ways of the game commands.

        We shall be limited in our ability to read various DNA files because
        there are a number of files that only function for their respective Data Base files.
        Additionally, the small.txt or large.txt must have already been read.
        so that it can read the various DNA data in the range of the those
        files ('small.txt' or 'large.txt').

----------------------------------------------------------------------*/

#include "ourvector.h" // for the use of ourvector from the hader file (ourvector.h)
#include <fstream>     // We will need the libary to read from the file
#include <sstream>     // for stringstream()

using namespace std;

// Keep the data that is being read from the file into on
// stuct that has all the ourvactor init
struct DnaProfileApp
{
    ourvector<string> dbData; // store the first line from the data base file

    ourvector<ourvector<char>> str; // store the Data base file data in dbData vac

    ourvector<char> dnaData; // store the DNA file data in dnaData vac

    ourvector<string> process; // store the index of the sequence of the after processing it

    ourvector<string> namesFromDb; // store the names of the people from the Data base file

    ourvector<int> foundIdx; // this will store the index of the number of the matching sequence
};

// Reading the Data base from the files
void loadDb(DnaProfileApp &dbFileData, string fileName)
{
    // read into line
    string line, names;
    ifstream readFile(fileName);

    // check file is open correctly
    if (readFile.fail())
    {
        cout << "Error: unable to open '" << fileName << "'" << endl;
        return;
    }

    getline(readFile, line, ','); // git rid of 'name'

    // get the STRs from the first line
    readFile >> line;

    ourvector<char> storeStrTmp; // this will store the first line string into char
                                 // and letter will be push it back in to str 2D vac

    // split the STR line and store it into Str Vac file
    for (unsigned i = 0; i < line.length(); i++)
    {
        if (line[i] == ',')
        {
            dbFileData.str.push_back(storeStrTmp);
            storeStrTmp.clear();
        }
        else
            storeStrTmp.push_back(line[i]);
    }

    dbFileData.str.push_back(storeStrTmp); // Push the last STR in vac

    // continue reading the file and push the content
    // from the file that is left reading
    while (readFile >> line)
    {
        stringstream ss(line);
        getline(ss, line, '\n'); // read up the new line
        dbFileData.dbData.push_back(line);
    }

    readFile.close();
}

// Load the DANs file that push that data into dna vac
void loadDna(DnaProfileApp &dnaFileData, string fileName)
{
    string line; // read into line

    ifstream readFile(fileName);

    // check if the file it valid to open
    if (readFile.fail())
    {
        cout << "Error: unable to open '" << fileName << "'" << endl;
        return;
    }
    // read the file and push the data into dna vac
    while (readFile >> line)
    {
        // split the line and push the data into DNA vac as a char
        for (unsigned i = 0; i < line.length(); i++)
        {
            dnaFileData.dnaData.push_back(line[i]);
        }
    }
}

// find the longes consecutive into the DNA data ourvactor
void findLongCons(DnaProfileApp &longCons, int start)
{
    // idx will be the first char from the str 2D vac that is
    // the going to be search in DNAs data
    // found will be incremented when the all the char matches with in the DANs data file
    // max will be the longes consecutive that will be push into 'process' vac
    int idx = 0, found = 0, max = 0;

    int tmpStart = 0; // tmpStart will have the last index of the char founded in the DNA data

    for (int finLCon = 0; finLCon < longCons.dnaData.size(); finLCon++)
    {
        // search the Char into DNAs data by comparing the current char of the DNA data
        // starting index and the first char of that starting index (idx)
        if (longCons.dnaData[finLCon] == longCons.str[start][idx])
        {
            // if idx is 0, meaning that the value is either reset it in between or the start of the search
            if (idx == 0)
            {
                tmpStart = finLCon;
            }
            idx++; // if char matches with str char
        }
        else
        {
            if (idx != 0)
            {
                finLCon = tmpStart;
            }
            // reset the value of the idx, and found. start idx from 0 so
            // that we start searching for the fist char into the DNA data reset the found
            idx = 0;
            found = 0;
        }

        if (idx == longCons.str[start].size()) // if the whole STR is found meaning that the all the char
                                               // of the str matches then we found the STR in the DNA data
        {
            // found goes up and idex goes back to 0
            found++;
            idx = 0;
            // check the found is more then current max if it is then the longes
            // consecutive will the the current founded
            if (found >= max)
            {
                max = found;
                tmpStart = finLCon;
            }
        }
    }
    longCons.process.push_back(to_string(max));
}

// When the input is process it goes into this function
// that checks validation first then move to counting the longes consecutive
void process(DnaProfileApp &process)
{

    // check before counting the longes consecutive
    if (process.dbData.size() == 0)
    {
        cout << "No database loaded.\n"
             << endl;
    }
    else if (process.dnaData.size() == 0)
    {
        cout << "No DNA loaded.\n"
             << endl;
    }

    // depending of the STRs ourvactor size it will call the
    // findLongCons() that many times and start at the first
    // till the size
    else
    {
        cout << "Processing DNA...\n";

        for (int i = 0; i < process.str.size(); i++)
        {
            findLongCons(process, i);
        }
    }
}

// This will search the longes consecutive that whom then belong to
void searchInDbForPerson(DnaProfileApp &searchInDb)
{
    string stringToBeSearch; // this will the string to be search

    // add comma to the string so that way we can easily seach the sequences count
    for (int i = 0; i < searchInDb.process.size(); i++)
    {
        string storeCurrCon = searchInDb.process[i];
        stringToBeSearch.append(storeCurrCon + ',');
    }
    // remove the names from the DB data vac and find the whom the longes sequences matches to
    for (int j = 0; j < searchInDb.dbData.size(); j++)
    {
        stringstream ss(searchInDb.dbData[j]);

        getline(ss, searchInDb.dbData[j], ',');
        searchInDb.namesFromDb.push_back(searchInDb.dbData[j]); // push the names into the names vac

        getline(ss, searchInDb.dbData[j], '\n');

        string currDbDataNumSeq = searchInDb.dbData[j] + ',';

        // check the the string we founded using the process function if
        // that string match with the current string of the DB data
        if (currDbDataNumSeq == stringToBeSearch)
        {
            searchInDb.foundIdx.push_back(j);
        }
    }
}

// seach the DNA sequence whom it belong to
void search(DnaProfileApp searchInDb)
{
    // will have the index that we founded if there are many then it will give us the
    // persons name whom it belongs to
    // check the validation first then if it is valid then move the else
    if (searchInDb.dbData.size() == 0)
        cout << "No database loaded.\n"
             << endl;

    else if (searchInDb.dnaData.size() == 0)
        cout << "No DNA loaded.\n"
             << endl;

    else if (searchInDb.process.size() == 0)
        cout << "No DNA processed.\n"
             << endl;

    else
    {
        cout << "Searching database...\n";
        // seach it if valid
        searchInDbForPerson(searchInDb);

        // if not found then
        if (searchInDb.foundIdx.size() == 0)
        {
            cout << "Not found in database.\n"
                 << endl;
        }
        // if found then print it out
        // if there are many matches then print out those names as well
        else
        {
            cout << "Found in database! DNA matches: ";
            for (int i = 0; i < searchInDb.foundIdx.size(); i++)
            {
                string name = searchInDb.namesFromDb[searchInDb.foundIdx[i]];
                cout << name << endl;
            }
        }
    }
}

// Display the contented of the Data base vac
void displayDB(DnaProfileApp &display)
{
    cout << "Database loaded: " << endl;

    // print the data
    for (auto i = 0; i < display.dbData.size(); i++)
    {
        for (unsigned j = 0; j < display.dbData[i].size(); j++)
        {
            if (display.dbData[i][j] != ',')
                cout << display.dbData[i][j];
            else
                cout << " ";
        }

        cout << endl;
    }
}

// display the longes convective that we have found it using process function
void displayLongesCon(DnaProfileApp &display)
{
    // displaying the longes consecutive to it appropriate STR
    // and if the comma is found then put ':' and print it
    cout << "\nDNA processed, STR counts: " << endl;
    for (int i = 0; i < display.str.size(); i++)
    {
        for (int j = 0; j < display.str[i].size(); j++)
        {
            cout << display.str[i][j];
        }
        cout << ": " << display.process[i] << endl;
    }
    cout << endl;
}

// Read the file that user entered in the range start to end and push that data
// into DNA data ourvac
void addMultiFiles(DnaProfileApp &multiDans, int &startFile, int &stopFile)
{
    string file = "";
    for (int i = startFile; i < stopFile; i++) // stop at the stopFile
    {
        file = to_string(i) + ".txt"; // convert the int into string
        loadDna(multiDans, file);     // and read the file
        file.clear();
    }
}

// Creative component and the multiple files to a DNA data ourvac
void loadMultiDANs(DnaProfileApp &multiDans, string fileName)
{
    int startFile = 0, stopFile = 0;
    // check the validations
    if (fileName == " ")
    {
        cout << "No Data Base Loaded." << endl;
        return;
    }
    else
    {
        // get the start file input and stopFile so we can read the
        // files in the range of those files

        cout << "Enter a Start File Number (i.e '1'): ";
        cin >> startFile;

        cout << "Enter a Stop File (i.e '5'): ";
        cin >> stopFile;

        // if small.txt then check the range of the file
        if (fileName == "small.txt" && (startFile >= 1 && stopFile <= 5))
        {
            addMultiFiles(multiDans, startFile, stopFile);
        }

        // if large.txt then check the range of the file
        else if (fileName == "large.txt" && (startFile >= 6 && stopFile <= 20))
        {
            addMultiFiles(multiDans, startFile, stopFile);
        }
        // if not out of range
        else
        {
            cout << "Out of Range." << endl;
        }
    }
}

// display all the things when the display command is given as input
void display(DnaProfileApp &display)
{
    // check all the validation before printing any thing
    // if anything is not present the give an appropriate message
    if (display.dbData.size() == 0)
        cout << "No database loaded." << endl;
    else // go into function and print it

        displayDB(display);
    if (display.dnaData.size() == 0)
        cout << "\nNo DNA loaded.\n"
             << endl;
    else
    {
        cout << "DNA loaded: " << endl;
        for (int j = 0; j < display.dnaData.size(); j++)
            cout << display.dnaData[j];
        cout << endl;
    }
    if (display.process.size() == 0)
        cout << "No DNA has been processed." << endl;
    else // go into function and print it
        displayLongesCon(display);
}

// runs all the commands here
int main()
{
    system("clear");
    cout << "Welcome to the DNA Profiling Application.\n";
    // get the command input and file input
    string commandEntry, fileEntryForDb = " ", fileEntryForDna;
    DnaProfileApp dnaProApp; // use the ourvactor that are in 'DnaProfileApp' struct to store the things
    cout << "Enter command or # to exit: ";
    while (commandEntry != "#")
    {
        cin >> commandEntry; // command entry
        // if the user wants to read the files
        if (commandEntry == "load_db")
        {
            cin >> fileEntryForDb;
            cout << "Loading database...\n";
            dnaProApp.dbData.clear();
            loadDb(dnaProApp, fileEntryForDb); // lode the DB file
        }
        // if the user wants to read the DNA files
        if (commandEntry == "load_dna")
        {
            cin >> fileEntryForDna;
            cout << "Loading DNA...\n";
            dnaProApp.dnaData.clear();
            loadDna(dnaProApp, fileEntryForDna); // lode the DNA file
        }
        if (commandEntry == "process")
            process(dnaProApp);
        if (commandEntry == "search")
            search(dnaProApp);

        if (commandEntry == "display")
            display(dnaProApp);

        if (commandEntry == "load_multi_dnas")
            loadMultiDANs(dnaProApp, fileEntryForDb);
        cout << "Enter command or # to exit: ";
    }
    return 0;
}
