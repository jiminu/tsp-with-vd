#ifndef _HEADEROFVDSFILE_H
#define _HEADEROFVDSFILE_H

namespace V {

namespace GeometryTier {


class HeaderOfVDSFile
{
private:
    rg_FLOAT m_version;
    char*    m_pdbFilename;
    char*    m_date;
    char*    m_author;
    char*    m_organization;
    char*    m_operatingSystem;
    rg_INT   m_computingAlgorithm;  //  1: edge_tracing,  2: region_expansion
    char*    m_OS;
    char*    m_CPU;
    char*    m_RAM;
    rg_INT   m_computingTime;

public:
    HeaderOfVDSFile();
    HeaderOfVDSFile( rg_FLOAT version,
                     char*    pdbFilename,
                     char*    date,
                     char*    author,
                     char*    organization,
                     char*    operatingSystem,
                     rg_INT   computingAlgorithm,
                     char*    OS,
                     char*    CPU,
                     char*    RAM,
                     rg_INT   computingTime );
    HeaderOfVDSFile(const HeaderOfVDSFile& header);
    ~HeaderOfVDSFile();


};

} // namespace GeometryTier

} // namespace V


#endif

