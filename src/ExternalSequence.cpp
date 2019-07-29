#include "ExternalSequence.h"

#include <stdio.h>		// sscanf
#include <cstring>		// strlen etc
#include <iomanip>		// std::setw etc

#include <algorithm>	// for std::max_element

#include <math.h>		// fabs etc

ExternalSequence::PrintFunPtr ExternalSequence::print_fun = &ExternalSequence::defaultPrint;
const int ExternalSequence::MAX_LINE_SIZE = 256;
const char ExternalSequence::COMMENT_CHAR = '#';

/***********************************************************/
ExternalSequence::ExternalSequence()
{
	m_blocks.clear();	// probably not needed.
	version_major=0;
	version_minor=0;
	version_revision=0;
	version_combined=0;
}


/***********************************************************/
ExternalSequence::~ExternalSequence(){}

/***********************************************************/
void ExternalSequence::print_msg(MessageType level, std::ostream& ss) {
	if (MSG_LEVEL>=level) {
#if defined(VXWORKS) || defined (BUILD_PLATFORM_LINUX)
		// we skip messages on the scanner platforms due to performanc limitations
		// we could trivially use UTRACE on newer scanners, but it is not compatible with older platforms
#else		
		std::ostringstream oss;
		oss.width(2*(level-1)); oss << "";
		oss << static_cast<std::ostringstream&>(ss).str();
		print_fun(oss.str().c_str());
#endif
	}
}

/***********************************************************/
bool ExternalSequence::load(std::string path)
{
	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "Reading external sequence files");

	char buffer[MAX_LINE_SIZE];
	char tmpStr[MAX_LINE_SIZE];

	// **********************************************************************************************************************
	// ************************ READ SHAPES ***********************************

	// Try single file mode (everything in .seq file)
	std::ifstream data_file;
	bool isSingleFileMode = true;
	std::string filepath = path;
	if (filepath.substr(filepath.size()-4) != std::string(".seq")) {
		filepath = path + PATH_SEPARATOR + "external.seq";
	}
	// Open in binary mode to ensure all end-of-line characters are processed
	data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
	data_file.seekg(0, std::ios::beg);

	if (!data_file.good())
	{
		// Try separate file mode (blocks.seq, events.seq, shapes.seq)
		filepath = path + PATH_SEPARATOR + "shapes.seq";
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);
		isSingleFileMode = false;
	}
	if (!data_file.good())
	{
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
		return false;
	}
	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Building index" );

	// Save locations of section tags
	buildFileIndex(data_file);

	// Read version section
	if (m_fileIndex.find("[VERSION]") != m_fileIndex.end()) {
		print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "decoding VERSION section");
		// Version is a recommended but not a compulsory section
		// very basic reading code, repeated keywords will overwrite previous values, no serious error checking
		data_file.seekg(m_fileIndex["[VERSION]"], std::ios::beg);
		skipComments(data_file,buffer);			// load up some data and ignore comments & empty lines
		while (data_file.good() && buffer[0]!='[')
		{
			//print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "buffer: \n" << buffer << std::endl );
			if (0==strncmp(buffer,"major",5)) {
				    if (1!=sscanf(buffer+5, "%d", &version_major)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode version_major");
					return false;
				}
			    print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "major=" << version_major);		
			} else if (0==strncmp(buffer,"minor",5)) {
				if (1!=sscanf(buffer+5, "%d", &version_minor)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode version_minor");
					return false;
				}
				print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "minor=" << version_minor);
			}
			else if (0==strncmp(buffer,"revision",8)) {
				if (1!=sscanf(buffer+8, "%d", &version_revision)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode version_revision \n" << buffer << std::endl );
					return false;
				}
				print_msg(DEBUG_MEDIUM_LEVEL, std::ostringstream().flush() << "revision=" << version_revision);
			}
			else
			{
				print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unknown field in the [VERSION] block");
				return false;
			}
			//getline(data_file, buffer, MAX_LINE_SIZE);
			skipComments(data_file,buffer);			// load up some data and ignore comments & empty lines
		}
		version_combined=version_major*1000000L+version_minor*1000L+version_revision;
	}
	
	// Read shapes section
	// ------------------------
	if (m_fileIndex.find("[SHAPES]") == m_fileIndex.end()) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: [SHAPES] section");
		return false;
	}
	data_file.seekg(m_fileIndex["[SHAPES]"], std::ios::beg);
	skipComments(data_file,buffer);			// Ignore comments & empty lines

	int shapeId, numSamples;
	float sample;

	m_shapeLibrary.clear();
	while (data_file.good() && buffer[0]=='s')
	{
		if (2!=sscanf(buffer, "%s%d", tmpStr, &shapeId)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode 'shapeId'\n" << buffer << std::endl );
			return false;
		}
		getline(data_file, buffer, MAX_LINE_SIZE);
		if (2!=sscanf(buffer, "%s%d", tmpStr, &numSamples)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode 'numSamples'\n" << buffer << std::endl );
			return false;
		}

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Reading shape " << shapeId );

		CompressedShape shape;
		shape.samples.clear();
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='s' || strlen(buffer)==0) {
				break;
			}
			if (1!=sscanf(buffer, "%f", &sample)) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode 'sample'\n" << buffer << std::endl );
				return false;
			}
			shape.samples.push_back(sample);
		}
		shape.numUncompressedSamples=numSamples;

		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape index " << shapeId << " has " << shape.samples.size()
			<< " compressed and " << shape.numUncompressedSamples << " uncompressed samples" );

		m_shapeLibrary[shapeId] = shape;

		skipComments(data_file,buffer);			// Ignore comments & empty lines
	}
	data_file.clear();	// In case EOF reached

	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- SHAPES READ numShapes: " << m_shapeLibrary.size() );


	// **********************************************************************************************************************
	// ************************ READ EVENTS ********************
	if (!isSingleFileMode) {
		filepath = path + PATH_SEPARATOR + "events.seq";
		data_file.close();
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);

		if (!data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}
		buildFileIndex(data_file);
	}

	// Read RF section
	// ------------------------
	if (m_fileIndex.find("[RF]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[RF]"], std::ios::beg);

		int rfId;
		m_rfLibrary.clear();
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			RFEvent event;
			if (version_combined<1002000L)
			{
				if (6!=sscanf(buffer, "%d%f%d%d%f%f", &rfId, &(event.amplitude),
							&(event.magShape),&(event.phaseShape),
							&(event.freqOffset), &(event.phaseOffset)
							)) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode RF event\n" << buffer << std::endl );
				return false;
			}
				event.delay=0;
			}
			else
			{
				if (7!=sscanf(buffer, "%d%f%d%d%d%f%f", &rfId, &(event.amplitude),
							&(event.magShape),&(event.phaseShape), &(event.delay),
							&(event.freqOffset), &(event.phaseOffset)
							)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode RF event\n" << buffer << std::endl );
					return false;
				}
			}
			m_rfLibrary[rfId] = event;
		}
	}
	
	// Read *arbitrary* gradient section
	// -------------------------------
	m_gradLibrary.clear();
	if (m_fileIndex.find("[GRADIENTS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[GRADIENTS]"], std::ios::beg);

		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			int gradId;
			GradEvent event;
			if ( version_combined>=1001001L )
			{
				if (4!=sscanf(buffer, "%d%f%d%d", &gradId, &(event.amplitude), &(event.shape), &(event.delay))) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode gradient event\n" << buffer << std::endl );
					return false;
				}
			}
			else
			{
				if (3!=sscanf(buffer, "%d%f%d", &gradId, &(event.amplitude), &(event.shape))) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode gradient event\n" << buffer << std::endl );
					return false;
				}
				event.delay=0;
			}
			m_gradLibrary[gradId] = event;
		}
	}

	// Read *trapezoid* gradient section
	// -------------------------------
	if (m_fileIndex.find("[TRAP]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[TRAP]"], std::ios::beg);

		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			int gradId;
			GradEvent event;
			if ( version_combined>=1001001L )
			{
				if (6!=sscanf(buffer, "%d%f%ld%ld%ld%d", &gradId, &(event.amplitude),
					&(event.rampUpTime),&(event.flatTime),&(event.rampDownTime),&(event.delay))) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode trapezoid gradient entry" << buffer << std::endl );
					return false;
				}					
			}
			else
			{
				if (5!=sscanf(buffer, "%d%f%ld%ld%ld", &gradId, &(event.amplitude),
							&(event.rampUpTime),&(event.flatTime),&(event.rampDownTime))) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode trapezoid gradient entry" << buffer << std::endl );
					return false;
				}
				event.delay=0;
			}
			
			event.shape=0;
			m_gradLibrary[gradId] = event;
		}
	}

	// Sort gradients based on index
	// -----------------------------
	//std::sort(m_gradLibrary.begin(),m_gradLibrary.end(),compareGradEvents);


	// Read ADC section
	// -------------------------------
	if (m_fileIndex.find("[ADC]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[ADC]"], std::ios::beg);

		int adcId;
		m_adcLibrary.clear();
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			ADCEvent event;
			if (6!=sscanf(buffer, "%d%d%d%d%f%f", &adcId, &(event.numSamples),
						&(event.dwellTime),&(event.delay),&(event.freqOffset),&(event.phaseOffset))) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode ADC event\n" << buffer << std::endl );
				return false;
			}
			m_adcLibrary[adcId] = event;
		}
	}

	// Read delays section
	// -------------------------------
	if (m_fileIndex.find("[DELAYS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[DELAYS]"], std::ios::beg);

		int delayId;
		long delay;
		m_delayLibrary.clear();
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			if (2!=sscanf(buffer, "%d%ld", &delayId, &delay)) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode delay event\n" << buffer << std::endl );
				return false;
			}
			m_delayLibrary[delayId] = delay;
		}
	}

	// Read extensions section
	// -------------------------------
	m_extensionLibrary.clear();
	m_extensionNameIDs.clear();
	m_triggerLibrary.clear(); // clear also all known extension libraries
	std::map<std::string,int>::iterator itFI = m_fileIndex.find("[EXTENSIONS]");
	if ( itFI != m_fileIndex.end()) {
		data_file.seekg(itFI->second, std::ios::beg);
		std::set<int>::iterator itSFI = m_fileSections.find(itFI->second);
		if ( itSFI==m_fileSections.end() ||
			 (++itSFI)==m_fileSections.end() )
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed find the end of the section while reading EXTENSIONS");
			return false;
		}
		int sectionEnd = *itSFI;
		// we first read in the extension list
		int nID;
		int nExtensionID=EXT_LIST; // EXT_LIST means we are reading the extension list
		while ( data_file.tellg()<sectionEnd &&
			    getline(data_file, buffer, MAX_LINE_SIZE)) 
		{
			if (buffer[0]=='#' || buffer[0]=='[' || strlen(buffer)==0) {
				continue;
			}
			print_msg(ERROR_MSG, std::ostringstream().flush() << "input line: " << buffer);
			if (0==strncmp(buffer,"extension",9)) {
				// read new extension ID from the header
				char szStrID[MAX_LINE_SIZE];
				int nInternalID=0;
				int nKnownID=EXT_UNKNOWN;
				if (2!=sscanf(buffer, "extension %s %d", szStrID, &nInternalID)) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode extension header entry\n" << buffer << std::endl );
					return false;
				}
				// here is the list if extensions we currently recognize
				if (0==strcmp("TRIGGERS",szStrID))
					nKnownID=EXT_TRIGGER;
				else if (0==strcmp("ROTATIONS",szStrID))
					nKnownID=EXT_ROTATION;
				if (nKnownID!=EXT_UNKNOWN)
					m_extensionNameIDs[nInternalID]=std::make_pair<std::string,int>(szStrID,nKnownID);
				else {
					print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unknown extension ignored\n" << buffer << std::endl );
				}
				nExtensionID=nKnownID;
			}
			else
			{
				ExtensionListEntry extEntry;
				TriggerEvent trigger;
				RotationEvent rotation;
				switch (nExtensionID) {
					case EXT_LIST: 
						if (4!=sscanf(buffer, "%d%d%d%ld", &nID, &(extEntry.type), &(extEntry.ref), &(extEntry.next))) {
							print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode extension list entry\n" << buffer << std::endl );
							return false;
						}
                        print_msg(ERROR_MSG, std::ostringstream().flush() << "decoding extension list entry " << buffer);
						m_extensionLibrary[nID] = extEntry;
						print_msg(ERROR_MSG, std::ostringstream().flush() << "nID:" << nID << " type:" << extEntry.type << " ref" << extEntry.ref << " next:" << extEntry.next);
						break;
					case EXT_TRIGGER: 
						if (5!=sscanf(buffer, "%d%d%d%ld%ld", &nID, &(trigger.triggerType), &(trigger.triggerChannel), &(trigger.delay), &(trigger.duration))) {
							print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode trigger event\n" << buffer << std::endl );
							return false;
						}
						m_triggerLibrary[nID] = trigger;
						break;
					case EXT_ROTATION: 
						if (10!=sscanf(buffer, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &nID, 
									&rotation.rotMatrix[0], &rotation.rotMatrix[1], &rotation.rotMatrix[2],
									&rotation.rotMatrix[3], &rotation.rotMatrix[4], &rotation.rotMatrix[5],
									&rotation.rotMatrix[6], &rotation.rotMatrix[7], &rotation.rotMatrix[8])) {
							print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode rotation event\n" << buffer << std::endl );
							return false;
						}
						rotation.defined=true;
						m_rotationLibrary[nID] = rotation; 
						break;
					case EXT_UNKNOWN:
						break; // just ignore unknown extensions
				}
			}
		}
	}
	
	// Read gradient rotation section
	// -------------------------------
	/*if (m_fileIndex.find("[ROTATIONS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[ROTATIONS]"], std::ios::beg);

		int controlId;
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			ControlEvent event;
			event.type = ControlEvent::ROTATION;
			if (10!=sscanf(buffer, "%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &controlId, 
						&event.rotMatrix[0], &event.rotMatrix[1], &event.rotMatrix[2],
						&event.rotMatrix[3], &event.rotMatrix[4], &event.rotMatrix[5],
						&event.rotMatrix[6], &event.rotMatrix[7], &event.rotMatrix[8])) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode rotation event\n" << buffer << std::endl );
				return false;
			}
			m_controlLibrary[controlId] = event;
		}
	}*/
	
	
	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- EVENTS READ: "
		<<" RF: " << m_rfLibrary.size()
		<<" GRAD: " << m_gradLibrary.size()
		<<" ADC: " << m_adcLibrary.size()
		<<" DELAY: " << m_delayLibrary.size()
		<<" EXTENSIONS: " << m_extensionLibrary.size() + m_triggerLibrary.size() + m_rotationLibrary.size() );


	// **********************************************************************************************************************
	// ************************ READ BLOCKS ********************
	if (!isSingleFileMode) {
		filepath = path + PATH_SEPARATOR + "blocks.seq";
		data_file.close();
		data_file.open(filepath.c_str(), std::ios::in | std::ios::binary);
		data_file.seekg(0, std::ios::beg);

		if (!data_file.good())
		{
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Failed to read file " << filepath);
			return false;
		}
		buildFileIndex(data_file);
	}


	// Read definition section
	// ------------------------
	unsigned int numBlocks = 0;
	if (m_fileIndex.find("[DEFINITIONS]") != m_fileIndex.end()) {
		data_file.seekg(m_fileIndex["[DEFINITIONS]"], std::ios::beg);

		// Read each definition line
		m_definitions.clear();
		while (getline(data_file, buffer, MAX_LINE_SIZE)) {
			if (buffer[0]=='[' || strlen(buffer)==0) {
				break;
			}
			std::istringstream ss(buffer);
			std::string key;
			ss >> key;
			double value;
			std::vector<double> values;
			while (ss >> value) {
				values.push_back(value);
			}
			m_definitions[key] = values;
		}

		std::ostringstream out;
		out << "-- " << "DEFINITIONS READ: " << m_definitions.size() << " : ";
		for (std::map<std::string,std::vector<double> >::iterator it=m_definitions.begin(); it!=m_definitions.end(); ++it)
		{
			out<< it->first << " ";
			for (int i=0; i<it->second.size(); i++)
				out << it->second[i] << " ";
		}

		print_msg(DEBUG_HIGH_LEVEL, out);

	} // if definitions exist

	// Read blocks section
	// ------------------------
	if (m_fileIndex.find("[BLOCKS]") == m_fileIndex.end()) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Required: [BLOCKS] section");
		return false;
	}
	data_file.seekg(m_fileIndex["[BLOCKS]"], std::ios::beg);

	int blockIdx;
	EventIDs events;

	// Read blocks
	m_blocks.clear();
	while (getline(data_file, buffer, MAX_LINE_SIZE)) {
		if (buffer[0]=='[' || strlen(buffer)==0) {
			break;
		}

		memset(events.id, 0, NUM_EVENTS*sizeof(int));

		int ret=sscanf(buffer, "%d%d%d%d%d%d%d%d", &blockIdx,
				&events.id[DELAY],                              // Delay
				&events.id[RF],                                 // RF
				&events.id[GX],&events.id[GY],&events.id[GZ],   // Gradients
				&events.id[ADC],                                // ADCs
				&events.id[EXT]                                 // Extensions
				);
		if (7>ret
				) {
					print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: failed to decode event table entry:\n" << buffer << std::endl );
					print_msg(ERROR_MSG, std::ostringstream().flush() << "***        number of fields read: " << ret << std::endl );
			return false;
		}

		if (!checkBlockReferences(events)) {
			print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Block " << blockIdx
				<< " contains references to undefined events" );
			print_msg(ERROR_MSG, std::ostringstream().flush() << "***        RF:" << events.id[RF] << " GX:" << events.id[GX] << " GY:" << events.id[GY] << " GZ:" << events.id[GZ] << " ADC:" << events.id[DELAY] << " GX:" << events.id[DELAY] << " EXT:" << events.id[EXT]);
			return false;
		}
		// Add event IDs to list of blocks
		m_blocks.push_back(events);
	}
	data_file.close();

	print_msg(DEBUG_HIGH_LEVEL, std::ostringstream().flush() << "-- BLOCKS READ: " << m_blocks.size());

	// Num_Blocks definition (if defined) is used to check the correct number of blocks are read
	if (numBlocks>0 && m_blocks.size()!=numBlocks) {
		print_msg(ERROR_MSG, std::ostringstream().flush() << "*** ERROR: Expected " << numBlocks
		    << " blocks but read " << m_blocks.size() << " blocks");
		return false;
	}
	std::vector<double> def = GetDefinition("Scan_ID");
	int scanID = def.empty() ? 0: (int)def[0];
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "===== EXTERNAL SEQUENCE #" << std::setw(5) << scanID << " ===========" );
	print_msg(NORMAL_MSG, std::ostringstream().flush() << "==========================================" );

	return true;
};


/***********************************************************/
void ExternalSequence::skipComments(std::ifstream &fileStream, char *buffer)
{
	while (getline(fileStream, buffer, MAX_LINE_SIZE)) {
		if (buffer[0]!=COMMENT_CHAR && strlen(buffer)>0) {
			break;
		}
	}
};


/***********************************************************/
void ExternalSequence::buildFileIndex(std::ifstream &fileStream)
{
	char buffer[MAX_LINE_SIZE];
	
	while (getline(fileStream, buffer, MAX_LINE_SIZE)) {
		std::string line = std::string(buffer);
		if (line[0]=='[' && line[line.length()-1]==']') {
			m_fileIndex[line] = fileStream.tellg();
			m_fileSections.insert(fileStream.tellg());			
		}
	}
	m_fileSections.insert(fileStream.tellg()); // add the end-of-file (+1?)
	fileStream.clear();		// reset EOF flag
	fileStream.seekg(0, std::ios::beg);
};

/***********************************************************/
SeqBlock*	ExternalSequence::GetBlock(int index) {
	SeqBlock *block = new SeqBlock();

	// Copy event IDs
	EventIDs events = m_blocks[index];
	std::copy(events.id,events.id+NUM_EVENTS,&block->events[0]);

	//ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "GetBlock(" << index << ") : [ " << events.id[0] << " " << events.id[1] << " " << events.id[2] << " " << events.id[3] << " " << events.id[4] << " " << events.id[5] << " " << events.id[6] << "  ]");

	// Set some defaults
	block->index = index;
	block->delay = 0;
	block->adc.numSamples = 0;
	block->trigger.triggerType=0;
	block->trigger.triggerChannel=0;
	block->rotation.defined=false;

	// Set event structures (if applicable) so e.g. gradient type can be determined
	if (events.id[RF]>0)     block->rf      = m_rfLibrary[events.id[RF]];
	if (events.id[ADC]>0)    block->adc     = m_adcLibrary[events.id[ADC]];
	if (events.id[DELAY]>0)  block->delay   = m_delayLibrary[events.id[DELAY]];
	//if (events.id[CTRL]>0)   block->control = m_controlLibrary[events.id[CTRL]];
	for (unsigned int i=0; i<NUM_GRADS; i++)
		if (events.id[GX+i]>0) block->grad[i] = m_gradLibrary[events.id[GX+i]];
	// unpack (known) extension objects
	if (events.id[EXT]>0) {
		// oh yeah, the current data stuctures seem to be really ugly and slow...
		int nNextExtID=events.id[EXT];
		while (nNextExtID) {
			std::map<int,ExtensionListEntry>::iterator itEL = m_extensionLibrary.find(nNextExtID);
			if (itEL == m_extensionLibrary.end()) {
				print_msg(ERROR_MSG, std::ostringstream().flush() << "ERROR: could not find extension list entry " << nNextExtID);
				//return NULL;
				break;
			}
			// attempt to recognize the extension reference
			std::map<int,std::pair<std::string,int> >::iterator itEN=m_extensionNameIDs.find(itEL->second.type);
			if (itEN!=m_extensionNameIDs.end()) {
				// we have a known extension
				switch (itEN->second.second) {
					case EXT_TRIGGER:
						if (block->trigger.triggerType!=0) {
							print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: only one trigger per block is supported; error block: " << index );
						}
						else {
							// ok, lets find the trigger in the library
							block->trigger=m_triggerLibrary[itEL->second.ref]; // do we have to check whether it can be found?
						}
						break;
					case EXT_ROTATION:
						if (block->rotation.defined) {
							print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: only one rotation per block is supported; error block: " << index );
						}
						else {
							// ok, lets find the trigger in the library
							block->rotation=m_rotationLibrary[itEL->second.ref]; // do we have to check whether it can be found?
						}
						break;
					default:
						print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unimplemented extension type " << itEN->second.first << " in block " << index );
				}
			}
			else
			{
				print_msg(WARNING_MSG, std::ostringstream().flush() << "*** WARNING: unrecognized extension type " << itEL->second.type << " in block " << index );
			}
			// update the next pointer
			nNextExtID=itEL->second.next;
		}
	}
	// Calculate duration of block
	long duration = 0;
	if (block->isRF()) {
		RFEvent &rf = block->GetRFEvent();
		duration = MAX(duration, rf.delay+(long)m_shapeLibrary[rf.magShape].numUncompressedSamples);
	}

	for (int iC=0; iC<NUM_GRADS; iC++)
	{
		GradEvent &grad = block->GetGradEvent(iC);
		if (block->isArbitraryGradient(iC))
			duration = MAX(duration, (long)(10*m_shapeLibrary[grad.shape].numUncompressedSamples) + grad.delay);
		else if (block->isTrapGradient(iC))
			duration = MAX(duration, grad.rampUpTime + grad.flatTime + grad.rampDownTime + grad.delay);
	}
	if (block->isADC()) {
		ADCEvent &adc = block->GetADCEvent();
		duration = MAX(duration, adc.delay + (adc.numSamples*adc.dwellTime)/1000);
	}
	if (block->isTrigger()) {
		TriggerEvent &trigger = block->GetTriggerEvent();
		duration = MAX(duration, trigger.delay+trigger.duration );
	}

	// handling of delays has changed in revision 1.2.0
	if (version_combined<1002000L)
		block->duration = duration + block->delay;
	else
		block->duration = MAX(duration, block->delay);

	//ExternalSequence::print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "block duration: " << block->duration);
    
	return block;
}

/***********************************************************/
bool ExternalSequence::decodeBlock(SeqBlock *block)
{
	int *events = &block->events[0];
	print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Decoding block " << block->index << " events: "
		<< events[0]+1 << " " << events[1]+1 << " " << events[2]+1 << " "
		<< events[3]+1 << " " << events[4]+1 );

	std::vector<float> waveform;

	// Decode RF
	if (block->isRF())
	{
		// Decompress the shape for this channel
		CompressedShape& shape = m_shapeLibrary[block->rf.magShape];
		waveform.resize(shape.numUncompressedSamples);
		if (!decompressShape(shape,&waveform[0]))
			return false;

		/*
		// detect zero-padding 
		int nStart;
		for (nStart=0; nStart<shape.numUncompressedSamples; ++nStart)
			if (0.0 != waveform[nStart])
				break;
		int nEnd;
		for (nEnd=shape.numUncompressedSamples; nEnd>0; --nEnd)
		{
			//print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Tail 0 detection: nEnd=" << nEnd << " waveform[nEnd-1]=" << waveform[nEnd-1] );
			if (fabs(waveform[nEnd-1]) > 1e-6)
				break;
		}
		// mz: now only copy the inner non-zero part
		//block->rfAmplitude = std::vector<float>(waveform);
		block->rfAmplitude = std::vector<float>(waveform.begin()+nStart, waveform.begin()+nEnd);
		// mz: and take notice of optimizations
		block->rf.rfZeroHead = nStart;
		block->rf.rfZeroTail = shape.numUncompressedSamples - nEnd;
		print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "shape length was "<< shape.numUncompressedSamples << "0 detection foung " << block->rf.rfZeroHead << " head zeros and " << block->rf.rfZeroTail << " trailing zeros, new length is " << block->rfAmplitude.size() );

		CompressedShape& shapePhase = m_shapeLibrary[block->rf.phaseShape];
		waveform.resize(shapePhase.numUncompressedSamples);
		decompressShape(shapePhase,&waveform[0]);

		// Scale phase by 2pi
		std::transform(waveform.begin(), waveform.end(), waveform.begin(), std::bind1st(std::multiplies<float>(),TWO_PI));
		//block->rfPhase = std::vector<float>(waveform);
		// mz: now only copy the inner non-zero part
		block->rfPhase = std::vector<float>(waveform.begin()+nStart, waveform.begin()+nEnd);
		*/

		//MZ: original Kelvin's code follows
		block->rfAmplitude = std::vector<float>(waveform);

		CompressedShape& shapePhase = m_shapeLibrary[block->rf.phaseShape];
		waveform.resize(shapePhase.numUncompressedSamples);
		if (!decompressShape(shapePhase,&waveform[0]))
			return false;

		// Scale phase by 2pi
		std::transform(waveform.begin(), waveform.end(), waveform.begin(), std::bind1st(std::multiplies<float>(),TWO_PI));
		block->rfPhase = std::vector<float>(waveform);
	}

	// Decode gradients
	for (int iC=GX; iC<ADC; iC++)
	{

		if (block->isArbitraryGradient(iC-GX))	// is arbitrary gradient?
		{
			// Decompress the arbitrary shape for this channel
			CompressedShape& shape = m_shapeLibrary[block->grad[iC-GX].shape];

			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Loaded shape with "
				<< shape.samples.size() << " compressed samples" );

			waveform.resize(shape.numUncompressedSamples);
			if (!decompressShape(shape,&waveform[0]))
				return false;

			print_msg(DEBUG_LOW_LEVEL, std::ostringstream().flush() << "Shape uncompressed to "
				<< shape.numUncompressedSamples << " samples" );

			block->gradWaveforms[iC-GX] = std::vector<float>(waveform);
		}
	}

	// Decode ADC
	/*if (block->isADC())
	{
		block->adc = m_adcLibrary[events[ADC]];
	}

	// Decode Delays
	if (block->isDelay())
	{
		block->delay = m_delayLibrary[events[DELAY]];
	}*/

	checkGradient(*block);
	checkRF(*block);

	return true;
}

/***********************************************************/
bool ExternalSequence::decompressShape(CompressedShape& encoded, float *shape)
{
	float *packed = &encoded.samples[0];
	int numPacked = encoded.samples.size();
	int numSamples = encoded.numUncompressedSamples;

	int countPack=1;
	int countUnpack=1;
	while (countPack<numPacked)
	{

		if (packed[countPack-1]!=packed[countPack])
		{
			shape[countUnpack-1]=((float)packed[countPack-1]);
			countPack++; countUnpack++;
		}
		else
		{
			int rep = ((int)packed[countPack+1])+2;
			if (fabs(packed[countPack+1]+2-rep)>1e-6) // MZ: detect format error present in some Pulseq Matlab toolbox versions
			{
				print_msg(ERROR_MSG, std::ostringstream().flush() << "ERROR: compressed shape format error detected \n"
																	 "  packed[countPack-1]=" << packed[countPack-1] << "  packed[countPack]=" << packed[countPack] << std::endl <<
																	 "  packed[countPack+1]=" << packed[countPack+1] << "  rep=" << rep << "  countPack=" << countPack );
				return false;
			}
			for (int i=countUnpack-1; i<=countUnpack+rep-2; i++)
				shape[i]= ((float)packed[countPack-1]);
			countPack += 3;
			countUnpack += rep;
		}

	}
	if (countPack==numPacked) {
		shape[countUnpack-1]=((float)packed[countPack-1]);
	}

	// Cumulative sum
	for (int i=1; i<numSamples; i++) {
		shape[i]=shape[i]+shape[i-1];
	}

	return true;
};

/***********************************************************/
bool ExternalSequence::checkBlockReferences(EventIDs& events)
{
	bool error;
	error = (events.id[RF]>0    && m_rfLibrary.count(events.id[RF])==0);
	error|= (events.id[GX]>0    && m_gradLibrary.count(events.id[GX])==0);
	error|= (events.id[GY]>0    && m_gradLibrary.count(events.id[GY])==0);
	error|= (events.id[GZ]>0    && m_gradLibrary.count(events.id[GZ])==0);
	error|= (events.id[ADC]>0   && m_adcLibrary.count(events.id[ADC])==0);
	error|= (events.id[DELAY]>0 && m_delayLibrary.count(events.id[DELAY])==0);
	//error|= (events.id[CTRL]>0  && m_controlLibrary.count(events.id[CTRL])==0); // TODO: currently all error checking is done in getBlock(); it needs to be done here 
	
	return (!error);
}

/***********************************************************/
void ExternalSequence::checkGradient(SeqBlock& block)
{
	for (unsigned int i=0; i<NUM_GRADS; i++)
	{
		std::vector<float> &waveform = block.gradWaveforms[i];
		for (unsigned j=0; j<waveform.size(); j++)
		{
			if (waveform[j]>1.0)  waveform[j]= 1.0;
			if (waveform[j]<-1.0) waveform[j]=-1.0;
		}
		// Ensure last point is zero // MZ: no, its wrong! trapezoid gradients have a non-zero at the end!
		// if (waveform.size()>0) waveform[waveform.size()-1]=0.0;
	}
}


/***********************************************************/
void ExternalSequence::checkRF(SeqBlock& block)
{
	for (unsigned int i=0; i<block.rfAmplitude.size(); i++)
	{
		if (block.rfAmplitude[i]>1.0) block.rfAmplitude[i]=1.0;
		if (block.rfAmplitude[i]<0.0) block.rfAmplitude[i]=0.0;
		if (block.rfPhase[i]>TWO_PI-1.e-4) block.rfPhase[i]=(float)(TWO_PI-1.e-4);
		if (block.rfPhase[i]<0) block.rfPhase[i]=0.0;
	}
}


/***********************************************************/
bool ExternalSequence::getline(std::istream& is, char *buffer, int MAX_SIZE)
{
	//std::cout << "ExternalSequence::getline()" << std::endl;

	int i=0;	// Current position in buffer
	for(;;) {
		char c = (char)is.get();

		switch (c) {
			case '\n':
				buffer[i]='\0';
				//std::cout << "ExternalSequence::getline() EOL i:" << i << std::endl;
				//std::cout << "buffer: " << buffer << std::endl;
				return true;
			case '\r':
				if(is.peek() == '\n') {
					is.get();       // Discard character
				}
				buffer[i]='\0';
				//std::cout << "ExternalSequence::getline() CR i:" << i << std::endl;
				//std::cout << "buffer: " << buffer << std::endl;
				return true;
			case EOF:
				if(i==0)
					is.seekg(-1);   // Create error on stream
				buffer[i]='\0';     // In case last line has no line ending
				//std::cout << "ExternalSequence::getline() EOF i:" << i << std::endl;
				return (i!=0);
			default:
				buffer[i++] = c;
				if (i>=MAX_SIZE)
				{
					return true;
				}
		}
	}
}
