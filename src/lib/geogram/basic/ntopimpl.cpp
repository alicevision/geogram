
#include <geogram/basic/logger.h>

namespace GEO {

Logger Logger::sLogger;

namespace Biblio {

/**
* \brief Initializes the bibliography system.
*/
void initialize() {}

/**
* \brief Terminates the bibliography system.
*/
void terminate() {}

/**
* \brief Registers a set of bibliographic references.
* \param[in] bib_refs a string with the bibliographic references,
*  in Bibtex format.
*/
void register_references(const char* bib_refs) {};

/**
* \brief Cites a bibliographic reference.
* \details Client code should not use this function and should use
*  the geo_cite() macro instead.
* \param[in] ref the citation key.
* \param[in] file the source filename from which the citation key is
*  cited.
* \param[in] line the source line number.
* \param[in] function the name of the function from which the citation
*  key is cited.
* \param[in] info more information about the context of the citation.
*/
void cite(
  const char* ref,
  const char* file, int line,
  const char* function,
  const char* info = nil
) {};

/**
* \brief Resets all citations.
*/
void reset_citations();
}

}
