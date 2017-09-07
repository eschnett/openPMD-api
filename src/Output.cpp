#include <ios>
#include <iostream>
#include <set>

#include "../include/Auxiliary.hpp"
#include "../include/Output.hpp"
#include "../include/IO/HDF5/HDF5IOHandler.hpp"
#include "../include/IO/HDF5/ParallelHDF5IOHandler.hpp"


char const * const Output::BASEPATH = "/data/%T/";
char const * const Output::OPENPMD = "1.0.1";

std::ostream&
operator<<(std::ostream& os, Output::IterationEncoding ie)
{
    switch( ie )
    {
        case Output::IterationEncoding::fileBased:
            os<<"fileBased";
            break;
        case Output::IterationEncoding::groupBased:
            os<<"groupBased";
            break;
    }
    return os;
}

Output::Output(std::string const& path,
               std::string const& name,
               IterationEncoding ie,
               Format f,
               AccessType at)
        : iterations{Container< Iteration, uint64_t >()},
          m_iterationEncoding{ie},
          m_name{name}
{
    switch( f )
    {
        case Format::HDF5:
            IOHandler = std::make_shared<HDF5IOHandler>(path, at);
            break;
        case Format::PARALLEL_HDF5:
            IOHandler = std::make_shared<ParallelHDF5IOHandler>(path, at);
            break;
        case Format::ADIOS:
            break;
        case Format::NONE:
            IOHandler = std::make_shared<NONEIOHandler>(path, at);
            break;
    }
    iterations.IOHandler = IOHandler;
    iterations.parent = this;
    switch( at )
    {
        case AccessType::CREAT:
        {
            setOpenPMD(OPENPMD);
            setOpenPMDextension(0);
            setAttribute("basePath", std::string(BASEPATH));
            setMeshesPath("meshes/");
            setParticlesPath("particles/");
            switch( ie )
            {
                case Output::IterationEncoding::fileBased:
                    setIterationFormat(m_name + "_%T");
                    setAttribute("iterationEncoding", std::string("fileBased"));
                    break;
                case Output::IterationEncoding::groupBased:
                    setIterationFormat("/data/%T/");
                    setAttribute("iterationEncoding", std::string("groupBased"));
                    break;
            }
            break;
        }
        case AccessType::READ_ONLY:
        {
            read();
            break;
        }
        case AccessType::READ_WRITE:
        {
            read();
            break;
        }
    }
}

Output::~Output()
{
    if( IOHandler->accessType == AccessType::CREAT || IOHandler->accessType == AccessType::READ_WRITE )
        flush();
    IOHandler->flush();
}

std::string
Output::openPMD() const
{
    return getAttribute("openPMD").get< std::string >();
}

Output&
Output::setOpenPMD(std::string const& o)
{
    setAttribute("openPMD", o);
    dirty = true;
    return *this;
}

uint32_t
Output::openPMDextension() const
{
    return getAttribute("openPMDextension").get< uint32_t >();
}

Output&
Output::setOpenPMDextension(uint32_t oe)
{
    setAttribute("openPMDextension", oe);
    dirty = true;
    return *this;
}

std::string
Output::basePath() const
{
    return getAttribute("basePath").get< std::string >();
}

Output&
Output::setBasePath(std::string const& bp)
{
    std::string version = openPMD();
    if( version == "1.0.0" || version == "1.0.1" )
        throw std::runtime_error("Custom basePath not available in openPMD <=1.0.1");

    setAttribute("basePath", bp);
    dirty = true;
    return *this;
}

std::string
Output::meshesPath() const
{
    return getAttribute("meshesPath").get< std::string >();
}

Output&
Output::setMeshesPath(std::string const& mp)
{
    if( ends_with(mp, "/") )
        setAttribute("meshesPath", mp);
    else
        setAttribute("meshesPath", mp + "/");
    dirty = true;
    return *this;
}

std::string
Output::particlesPath() const
{
    return getAttribute("particlesPath").get< std::string >();
}

Output&
Output::setParticlesPath(std::string const& pp)
{
    if( ends_with(pp, "/") )
        setAttribute("particlesPath", pp);
    else
        setAttribute("particlesPath", pp + "/");
    dirty = true;
    return *this;
}

std::string
Output::author() const
{
    return getAttribute("author").get< std::string >();
}

Output&
Output::setAuthor(std::string const& a)
{
    setAttribute("author", a);
    dirty = true;
    return *this;
}

std::string
Output::software() const
{
    return getAttribute("software").get< std::string >();
}

Output&
Output::setSoftware(std::string const& s)
{
    setAttribute("software", s);
    dirty = true;
    return *this;
}

std::string
Output::softwareVersion() const
{
    return getAttribute("softwareVersion").get< std::string >();
}

Output&
Output::setSoftwareVersion(std::string const& sv)
{
    setAttribute("softwareVersion", sv);
    dirty = true;
    return *this;
}

std::string
Output::date() const
{
    return getAttribute("date").get< std::string >();
}

Output&
Output::setDate(std::string const& d)
{
    setAttribute("date", d);
    dirty = true;
    return *this;
}

Output::IterationEncoding
Output::iterationEncoding() const
{
    return m_iterationEncoding;
}

Output&
Output::setIterationEncoding(IterationEncoding ie)
{
    if( written )
        throw std::runtime_error("A files iterationEncoding can not (yet) be changed after it has been written.");

    switch( ie )
    {
        case Output::IterationEncoding::fileBased:
            setIterationFormat(m_name + "_%T");
            setAttribute("iterationEncoding", std::string("fileBased"));
            break;
        case Output::IterationEncoding::groupBased:
            setIterationFormat("/data/%T/");
            setAttribute("iterationEncoding", std::string("groupBased"));
            break;
    }
    dirty = true;
    return *this;
}

std::string
Output::iterationFormat() const
{
    return getAttribute("iterationFormat").get< std::string >();
}

Output&
Output::setIterationFormat(std::string const& i)
{
    if( written )
        throw std::runtime_error("A files iterationFormat can not (yet) be changed after it has been written.");

    if( m_iterationEncoding == IterationEncoding::groupBased )
    {
        if( basePath() != i && (openPMD() == "1.0.1" || openPMD() == "1.0.0") )
            throw std::invalid_argument("iterationFormat must not differ from basePath " + basePath());
    }
    setAttribute("iterationFormat", i);
    dirty = true;
    return *this;
}

std::string
Output::name() const
{
    return m_name;
}

Output&
Output::setName(std::string const& n)
{
    if( written )
        throw std::runtime_error("A files name can not (yet) be changed after it has been written.");

    m_name = n;
    dirty = true;
    return *this;
}

void
Output::flush()
{
    if( IOHandler->accessType == AccessType::READ_WRITE || IOHandler->accessType == AccessType::CREAT )
    {
        switch( m_iterationEncoding )
        {
            using IE = IterationEncoding;
            case IE::fileBased:
            {
                for( auto& i : iterations )
                {
                    if( !i.second.written )
                        i.second.parent = this;
                    i.second.flushFileBased(i.first);

                    // TODO same problem as below
                    // [HDF5 backend: Container only corresponds to the ID written the earliest]
                    iterations.flush(replace_first(basePath(), "%T/", ""));

                    if( dirty )
                    {
                        Parameter< Operation::WRITE_ATT > attribute_parameter;
                        for( std::string const & att_name : attributes() )
                        {
                            attribute_parameter.name = att_name;
                            attribute_parameter.resource = getAttribute(att_name).getResource();
                            attribute_parameter.dtype = getAttribute(att_name).dtype;
                            // TODO "this" is too general for writing file-based root attributes
                            // [HDF5 backend: only one fileID is saved per Writable,
                            //  thus this Output only corresponds to the ID written the earliest]
                            IOHandler->enqueue(IOTask(this, attribute_parameter));
                        }
                        IOHandler->flush();
                        dirty = false;
                    }
                }
                break;
            }
            case IE::groupBased:
            {
                if( !written )
                {
                    Parameter< Operation::CREATE_FILE > file_parameter;
                    file_parameter.name = m_name;
                    IOHandler->enqueue(IOTask(this, file_parameter));
                    IOHandler->flush();
                }

                if( !iterations.written )
                    iterations.parent = this;
                iterations.flush(replace_first(basePath(), "%T/", ""));

                for( auto& i : iterations )
                {
                    if( !i.second.written )
                        i.second.parent = &iterations;
                    i.second.flushGroupBased(i.first);
                }

                flushAttributes();
            }
                break;
        }
    }
}

void
Output::read()
{
    //TODO find solution for handling fileBased output

    Parameter< Operation::OPEN_FILE > file_parameter;
    file_parameter.name = m_name;
    IOHandler->enqueue(IOTask(this, file_parameter));
    IOHandler->flush();
    /* allow all attributes to be set */
    written = false;

    using DT = Datatype;
    Parameter< Operation::READ_ATT > attribute_parameter;

    attribute_parameter.name = "openPMD";
    IOHandler->enqueue(IOTask(this, attribute_parameter));
    IOHandler->flush();
    if( *attribute_parameter.dtype == DT::STRING )
        setOpenPMD(strip(Attribute(*attribute_parameter.resource).get< std::string >(), {'\0'}));
    else
        throw std::runtime_error("Unexpected Attribute datatype for 'openPMD'");

    attribute_parameter.name = "openPMDextension";
    IOHandler->enqueue(IOTask(this, attribute_parameter));
    IOHandler->flush();
    if( *attribute_parameter.dtype == DT::UINT32 )
        setOpenPMDextension(Attribute(*attribute_parameter.resource).get< uint32_t >());
    else
        throw std::runtime_error("Unexpected Attribute datatype for 'openPMDextension'");

    attribute_parameter.name = "basePath";
    IOHandler->enqueue(IOTask(this, attribute_parameter));
    IOHandler->flush();
    if( *attribute_parameter.dtype == DT::STRING )
        setAttribute("basePath", strip(Attribute(*attribute_parameter.resource).get< std::string >(), {'\0'}));
    else
        throw std::runtime_error("Unexpected Attribute datatype for 'basePath'");

    attribute_parameter.name = "meshesPath";
    IOHandler->enqueue(IOTask(this, attribute_parameter));
    IOHandler->flush();
    if( *attribute_parameter.dtype == DT::STRING )
        setMeshesPath(strip(Attribute(*attribute_parameter.resource).get< std::string >(), {'\0'}));
    else
        throw std::runtime_error("Unexpected Attribute datatype for 'meshesPath'");

    attribute_parameter.name = "particlesPath";
    IOHandler->enqueue(IOTask(this, attribute_parameter));
    IOHandler->flush();
    if( *attribute_parameter.dtype == DT::STRING )
        setParticlesPath(strip(Attribute(*attribute_parameter.resource).get< std::string >(), {'\0'}));
    else
        throw std::runtime_error("Unexpected Attribute datatype for 'particlesPath'");

    attribute_parameter.name = "iterationEncoding";
    IOHandler->enqueue(IOTask(this, attribute_parameter));
    IOHandler->flush();
    if( *attribute_parameter.dtype == DT::STRING )
    {
        std::string encoding = strip(Attribute(*attribute_parameter.resource).get< std::string >(), {'\0'});
        if( encoding == "fileBased" )
        {
            if( m_iterationEncoding != IterationEncoding::fileBased)
                throw std::runtime_error("Supplied and read iterationEncoding do not match");
        }
        else if( encoding == "groupBased" )
        {
            if( m_iterationEncoding != IterationEncoding::groupBased)
                throw std::runtime_error("Supplied and read iterationEncoding do not match");
        } else
            throw std::runtime_error("Unknown iterationEncoding: " + encoding);
        setAttribute("iterationEncoding", encoding);
    }
    else
        throw std::runtime_error("Unexpected Attribute datatype for 'iterationEncoding'");

    attribute_parameter.name = "iterationFormat";
    IOHandler->enqueue(IOTask(this, attribute_parameter));
    IOHandler->flush();
    if( *attribute_parameter.dtype == DT::STRING )
        setIterationFormat(strip(Attribute(*attribute_parameter.resource).get< std::string >(), {'\0'}));
    else
        throw std::runtime_error("Unexpected Attribute datatype for 'iterationFormat'");

    iterations.clear();
    Parameter< Operation::OPEN_PATH > path_parameter;
    std::string version = openPMD();
    if( version == "1.0.0" || version == "1.0.1" )
        path_parameter.path = replace_first(basePath(), "/%T/", "");
    else
        throw std::runtime_error("Unknown openPMD version - " + version);
    IOHandler->enqueue(IOTask(&iterations, path_parameter));
    IOHandler->flush();

    iterations.readAttributes();

    /* obtain all paths inside the basepath (i.e. all iterations) */
    Parameter< Operation::LIST_PATHS > list_parameter;
    IOHandler->enqueue(IOTask(&iterations, list_parameter));
    IOHandler->flush();

    for( auto const& it : *list_parameter.paths )
    {
        Iteration& i = iterations[std::stoull(it)];
        path_parameter.path = it;
        IOHandler->enqueue(IOTask(&i, path_parameter));
        IOHandler->flush();
        i.read();
    }

    /* this file need not be flushed */
    iterations.written = true;
    written = true;
}