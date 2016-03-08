"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import rdflib

import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol


class AbstractPhenotypeAssociationSet(datamodel.DatamodelObject):
    compoundIdClass = datamodel.PhenotypeAssociationSetCompoundId

    def __init__(self, parentContainer, localId):
        super(AbstractPhenotypeAssociationSet, self).__init__(
            parentContainer, localId)

    def toProtocolElement(self):
        pas = protocol.PhenotypeAssociationSet()
        pas.name = self.getLocalId()
        pas.id = self.getId()
        pas.datasetId = self.getParentContainer().getId()
        pas.info = self.getInfo()
        return pas

    def getInfo(self):
        return {}


class SimulatedPhenotypeAssociationSet(AbstractPhenotypeAssociationSet):
    def __init__(self, parentContainer, localId, randomSeed):
        super(SimulatedPhenotypeAssociationSet, self).__init__(
            parentContainer, localId)

    def getAssociations(
            self, feature=None, environment=None, phenotype=None,
            pageSize=None, offset=0):
        if feature or environment or phenotype:
            fpa = protocol.FeaturePhenotypeAssociation()
            fpa.id = "test"
            fpa.features = []
            fpa.evidence = []
            fpa.environmentalContexts = []
            fpa.phenotype = protocol.PhenotypeInstance()
            fpa.phenotype.type = protocol.OntologyTerm()
            fpa.phenotype.type.id = "test"
            fpa.phenotypeAssociationSetId = self.getId()
            return [fpa]
        else:
            return []


class PhenotypeAssociationSet(AbstractPhenotypeAssociationSet):
    """
    An rdf object store.  The cancer genome database
    [Clinical Genomics Knowledge Base]
    (http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl),
    published by the Monarch project, was the source of Evidence.
    """

    def __init__(self, parentContainer, localId, dataDir):
        super(PhenotypeAssociationSet, self).__init__(parentContainer, localId)
        """
        Initialize dataset, using the passed dict of sources
        [{source,format}] see rdflib.parse() for more
        If path is set, this backend will load itself
        """

        # initialize graph
        self._rdfGraph = rdflib.ConjunctiveGraph()
        try:
            self._scanDataFiles(dataDir, ['*.ttl', '*.xml'])
        except AttributeError:
            pass

        # extract version
        cgdTTL = rdflib.URIRef("http://data.monarchinitiative.org/ttl/cgd.ttl")
        versionInfo = rdflib.URIRef(
            u'http://www.w3.org/2002/07/owl#versionInfo')
        self._version = None
        for s, p, o in self._rdfGraph.triples((cgdTTL, versionInfo, None)):
            self._version = o.toPython()

    def getAssociations(self, feature=None, environment=None,
                        phenotype=None, pageSize=None, offset=0):
        """
        This query is the main search mechanism.
        It queries the graph for annotations that match the
        AND of [feature,environment,phenotype].
        """
        # query to do search
        query = self._formatFilterQuery(feature, environment, phenotype)
        # DBG print(query)
        associations = self._rdfGraph.query(query)
        # DBG print(associations.bindings)
        # associations is now a dict with the following
# DBG TODO GS - please read this block and document as English comment
# [
#   {
#    rdflib.term.Variable(u'environment'):
#     rdflib.term.URIRef(u'http://www.drugbank.ca/drugs/DB00619'),
#    rdflib.term.Variable(u'feature'):
#     rdflib.term.URIRef(u'http://ohsu.edu/cgd/27d2169c'),
#    rdflib.term.Variable(u'environment_label'):
#     rdflib.term.Literal(u'imatinib'),
#    rdflib.term.Variable(u'feature_label'):
#     rdflib.term.Literal(u'KIT  wild type no mutation'),
#    rdflib.term.Variable(u'association'):
#     rdflib.term.URIRef(u'http://ohsu.edu/cgd/4657f28c'),
#    rdflib.term.Variable(u'phenotype_label'):
#     rdflib.term.Literal(u'GIST with decreased sensitivity to therapy'),
#    rdflib.term.Variable(u'sources'):
#     rdflib.term.Literal(u'http://www.ncbi.nlm.nih.gov/pubmed/18955458'),
#    rdflib.term.Variable(u'phenotype'):
#     rdflib.term.URIRef(u'http://ohsu.edu/cgd/87752f6c')
#   }
# ]

        # given get the details for the feature,phenotype and environment
        associations_details = self._detailTuples(
                                    self._extractAssociationsDetails(
                                        associations))
        # association_details is now a list of {subject,predicate,object}
        # for each of the association detail
# [ { u'object': u'http://purl.obolibrary.org/obo/SO_0001059',
#     u'predicate': u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type',
#     u'subject': u'http://ohsu.edu/cgd/27d2169c'},
#     .... ]

        HAS_QUALITY = 'http://purl.obolibrary.org/obo/BFO_0000159'
        associationList = []
        for association in associations.bindings:
            if '?feature' in association:
                association = self._bindingsToDict(association)
                association['feature'] = self._getDetails(
                                            association['feature'],
                                            associations_details)
                association['environment'] = self._getDetails(
                                        association['environment'],
                                        associations_details)
                association['phenotype'] = self._getDetails(
                                              association['phenotype'],
                                              associations_details)
                association['evidence'] = association['phenotype'][HAS_QUALITY]
                association['id'] = association['association']
                associationList.append(association)
# TODO - GS Please read and create English comment
# our association list is now a list of dicts
# [ {
# u'association':
#  u'http://ohsu.edu/cgd/4657f28c',
#     u'environment':
#      { u'http://purl.obolibrary.org/obo/RO_0002606':
#         u'http://ohsu.edu/cgd/b40a93d7',
#        u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
#         u'http://www.w3.org/2002/07/owl#Class',
#        u'http://www.w3.org/2000/01/rdf-schema#label':
#         u'imatinib',
#        u'http://www.w3.org/2000/01/rdf-schema#subClassOf':
#         u'http://purl.obolibrary.org/obo/CHEBI_23888',
#        u'id':
#         u'http://www.drugbank.ca/drugs/DB00619'
#      },
#     u'environment_label':
#       u'imatinib',
#     u'evidence':
#       u'http://ohsu.edu/cgd/decreased_sensitivity',
#     u'feature':
#       { u'http://purl.obolibrary.org/obo/RO_0002200':
#           u'http://ohsu.edu/cgd/87752f6c',
#         u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
#           u'http://purl.obolibrary.org/obo/SO_0001059',
#         u'http://www.w3.org/2000/01/rdf-schema#label':
#           u'KIT  wild type no mutation',
#         u'id':
#           u'http://ohsu.edu/cgd/27d2169c'
#       },
#     u'feature_label':
#       u'KIT  wild type no mutation',
#     u'phenotype':
#       { u'http://purl.obolibrary.org/obo/BFO_0000159':
#           u'http://ohsu.edu/cgd/decreased_sensitivity',
#         u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
#           u'http://purl.obolibrary.org/obo/OMIM_606764',
#         u'http://www.w3.org/2000/01/rdf-schema#label':
#           u'GIST with decreased sensitivity to therapy',
#         u'id':
#           u'http://ohsu.edu/cgd/87752f6c'},
#         u'phenotype_label':
#           u'GIST with decreased sensitivity to therapy',
#     u'sources':
#       u'http://www.ncbi.nlm.nih.gov/pubmed/18955458'
#  }]

        # create GA4GH objects
        associations = []
        for association in associationList:
            associations.append(self._toGA4GH(association))

        return associations

    def _extractAssociationsDetails(self, associations):
        """
        given a set of results from our search query, return a
        list of the URIRef for each `detail` (feature,environment,phenotype)
        """
        detailedURIRef = []
        for row in associations.bindings:
            # empty set [{}]
            if 'feature' in row:
                detailedURIRef.append(row['feature'])
                detailedURIRef.append(row['environment'])
                detailedURIRef.append(row['phenotype'])

# TESTING ... create test data
# obo = rdflib.Namespace('http://purl.obolibrary.org/obo/')
# oban = rdflib.Namespace('http://purl.org/oban/')
# dc = rdflib.Namespace('http://purl.org/dc/elements/1.1/')
# cgd = rdflib.Namespace('http://ohsu.edu/cgd/')
# faldo = rdflib.Namespace('http://biohackathon.org/resource/faldo#')
# drugBank = rdflib.Namespace('http://www.drugbank.ca/drugs/')
#
# namespace_manager = rdflib.namespace.NamespaceManager(rdflib.Graph())
# namespace_manager.bind('OBO', obo, override=False)
# namespace_manager.bind('OBAN', oban, override=False)
# namespace_manager.bind('CGD', cgd, override=False)
# namespace_manager.bind('dc', dc, override=False)
# namespace_manager.bind('faldo', faldo, override=False)
# namespace_manager.bind('DrugBank', drugBank, override=False)
#
# g = rdflib.Graph()
# g.namespace_manager = namespace_manager
# for row in associations:
#     for s, p, o in self._rdfGraph.triples((row['association'], None, None)):
#         g.add((s, p, o))
#
# for uriRef in detailedURIRef:
#     for s, p, o in self._rdfGraph.triples((uriRef, None, None)):
#         g.add((s, p, o))
# print(g.serialize(format='n3'))

        return detailedURIRef

    def _detailTuples(self, uriRefs):
        """
        given a list of uriRefs, return a list of dicts:
        {'subject': s, 'predicate': p, 'object': o }
        all values are strings
        """
        details = []
        for uriRef in uriRefs:
            for s, p, o in self._rdfGraph.triples((uriRef, None, None)):
                details.append({
                    'subject': s.toPython(),
                    'predicate': p.toPython(),
                    'object': o.toPython()
                })
        return details

    def _bindingsToDict(self, bindings):
        """
        given a binding from the sparql query result,
        create a dict of plain text
        """
        myDict = {}
        for k, v in bindings.iteritems():
            myDict[k.toPython().replace('?', '')] = v.toPython()
        return myDict

    def _addDataFile(self, filename):
        """ given a filename, add it to the graph """
        if filename.endswith('.ttl'):
            self._rdfGraph.parse(filename, format='n3')
        else:
            self._rdfGraph.parse(filename, format='xml')

    def _getDetails(self, uriRef, associations_details):
        """
        given a uriRef, return a dict of all the details for that Ref
        use the uriRef as the 'id' of the dict
        """
        associationDetail = {}
        for detail in associations_details:
            if detail['subject'] == uriRef:
                associationDetail[detail['predicate']] = detail['object']
            associationDetail['id'] = uriRef
        return associationDetail

    def _formatFilterQuery(self, feature=None, environment=None,
                           phenotype=None):
        """
        Generate a formatted sparql query with appropriate filters
        """
        query = """
            PREFIX OBAN: <http://purl.org/oban/>
            PREFIX OBO: <http://purl.obolibrary.org/obo/>
            PREFIX dc: <http://purl.org/dc/elements/1.1/>
            PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            SELECT
                ?association
                ?environment
                ?environment_label
                ?feature
                ?feature_label
                ?phenotype
                ?phenotype_label
                (GROUP_CONCAT(?source; separator="|") AS ?sources)
                ?evidence_type
                WHERE {
                    ?association  a OBAN:association .
                    ?association    OBO:RO_0002558 ?evidence_type .
                    ?association    OBO:RO_has_environment ?environment   .
                    OPTIONAL { ?association  dc:source ?source } .
                    ?association    OBAN:association_has_subject ?feature .
                    ?association    OBAN:association_has_object ?phenotype .
                    ?environment  rdfs:label ?environment_label  .
                    ?phenotype  rdfs:label ?phenotype_label  .
                    ?feature  rdfs:label ?feature_label  .
                    #%FILTER%
                    }
            GROUP  BY ?association
            ORDER  BY ?association
"""
        if feature is None and environment is None and phenotype is None:
            # TODO is this really the exception we want to throw?
            raise exceptions.NotImplementedException(
                "At least one of [feature, environment, phenotype] "
                "must be specified")
        filters = []

        # Strings
        if feature and isinstance(feature, basestring):
            filters.append('regex(?feature_label, "{}")'.format(feature))
        if environment and isinstance(environment, basestring):
            filters.append(
                'regex(?environment_label, "{}")'.format(environment))
        if phenotype and isinstance(phenotype, basestring):
            filters.append('regex(?phenotype_label, "{}")'.format(phenotype))

        # feature
        # ExternalIdentifier
        featureClause = ""
        if isinstance(feature, dict) and feature.get('ids'):
            features = []
            for _id in feature['ids']:
                features.append('?feature = <{}> '.format(
                    _id['database'] + _id['identifier']))
            featureClause = "({})".format(" || ".join(features))
            filters.append(featureClause)
        # OntologyTerms
        if isinstance(feature, dict) and feature.get('terms'):
            features = []
            for _term in feature['terms']:
                if _term.get('id'):
                    features.append('?feature = <{}> '.format(
                        _term['id']))
                else:
                    features.append('?feature = <{}> '.format(
                        self._toNamespaceURL(_term['term'])))
            featureClause = "({})".format(" || ".join(features))
            filters.append(featureClause)

        # environment
        # ExternalIdentifier
        if isinstance(environment, dict) and environment.get('ids'):
            environments = []
            for _id in environment['ids']:
                environments.append('?environment = <{}> '.format(
                    _id['database'] + _id['identifier']))
            environmentsClause = "({})".format(" || ".join(environments))
            filters.append(environmentsClause)

        # OntologyTerms
        if isinstance(environment, dict) and environment.get('terms'):
            environments = []
            for _term in environment['terms']:
                if _term.get('id'):
                    environments.append(
                        '?environment = <{}> '.format(_term['id']))
                else:
                    environments.append('?environment = <{}> '.format
                                        (self._toNamespaceURL(_term['term'])))
            environmentsClause = "({})".format(" || ".join(environments))
            filters.append(environmentsClause)

        # phenotype
        # ExternalIdentifier
        if isinstance(phenotype, dict) and phenotype.get('ids'):
            phenotypes = []
            for _id in phenotype['ids']:
                phenotypes.append('?phenotype = <{}> '.format(
                    _id['database'] + _id['identifier']))
            phenotypesClause = "({})".format(" || ".join(phenotypes))
            filters.append(phenotypesClause)
        # OntologyTerms
        if isinstance(phenotype, dict) and phenotype.get('terms'):
            phenotypes = []
            for _term in phenotype['terms']:
                if _term.get('id'):
                    phenotypes.append('?phenotype = <{}> '.format(
                        _term['id']))
                else:
                    phenotypes.append('?phenotype = <{}> '.format(
                        self._toNamespaceURL(_term['term'])))
            phenotypesClause = "({})".format(" || ".join(phenotypes))
            filters.append(phenotypesClause)

        # apply filters
        filter = "FILTER ({})".format(' && '.join(filters))
        query = query.replace("#%FILTER%", filter)
        return query

    def _toNamespaceURL(self, term):
        """
        Given an ontologyterm.term return namespace identifier.
        Leverages prefixes already in graph namespace
        Ex.  "DrugBank:DB01268" -> "http://www.drugbank.ca/drugs/DB01268"
        """
        (termPrefix, termId) = term.split(':')
        for prefix, namespace in self._rdfGraph.namespaces():
            if prefix == termPrefix:
                return namespace + termId
        raise exceptions.NotImplementedException(
           "Term has a prefix not found in this instance. {}".format(term))

    def _getIdentifier(self, url):
        """
        Given a url identifier return identifier portion
        Leverages prefixes already in graph namespace
        Returns None if no match
        Ex.  "http://www.drugbank.ca/drugs/DB01268" -> "DB01268"
        """
        for prefix, namespace in self._rdfGraph.namespaces():
            if namespace in url:
                return(url.replace(namespace, ''))

    def _getPrefix(self, url):
        """
        Given a url return namespace prefix.
        Leverages prefixes already in graph namespace
        Ex.  "http://www.drugbank.ca/drugs/" -> "Drugbank"
        """
        for prefix, namespace in self._rdfGraph.namespaces():
            if namespace.toPython() == url or namespace == url:
                return prefix
        raise exceptions.NotImplementedException(
           "No namespace found for url {}".format(url))

    def _getPrefixURL(self, url):
        """
        Given a url return namespace prefix.
        Leverages prefixes already in graph namespace
        Ex.  "http://www.drugbank.ca/drugs/DDD"
            -> "http://www.drugbank.ca/drugs/"
        """
        for prefix, namespace in self._rdfGraph.namespaces():
            if namespace.toPython() in url:
                return(namespace)

    def _toGA4GH(self, association):
        """
        given an association dict,
        return a protocol.FeaturePhenotypeAssociation
        """
        fpa = None
# DBG TODO GS - please read this block and document as English comment
# [ {
# u'association':
#  u'http://ohsu.edu/cgd/4657f28c',
#     u'environment':
#      { u'http://purl.obolibrary.org/obo/RO_0002606':
#         u'http://ohsu.edu/cgd/b40a93d7',
#        u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
#         u'http://www.w3.org/2002/07/owl#Class',
#        u'http://www.w3.org/2000/01/rdf-schema#label':
#         u'imatinib',
#        u'http://www.w3.org/2000/01/rdf-schema#subClassOf':
#         u'http://purl.obolibrary.org/obo/CHEBI_23888',
#        u'id':
#         u'http://www.drugbank.ca/drugs/DB00619'
#      },
#     u'environment_label':
#       u'imatinib',
#     u'evidence':
#       u'http://ohsu.edu/cgd/decreased_sensitivity',
#     u'feature':
#       { u'http://purl.obolibrary.org/obo/RO_0002200':
#           u'http://ohsu.edu/cgd/87752f6c',
#         u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
#           u'http://purl.obolibrary.org/obo/SO_0001059',
#         u'http://www.w3.org/2000/01/rdf-schema#label':
#           u'KIT  wild type no mutation',
#         u'id':
#           u'http://ohsu.edu/cgd/27d2169c'
#       },
#     u'feature_label':
#       u'KIT  wild type no mutation',
#     u'phenotype':
#       { u'http://purl.obolibrary.org/obo/BFO_0000159':
#           u'http://ohsu.edu/cgd/decreased_sensitivity',
#         u'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
#           u'http://purl.obolibrary.org/obo/OMIM_606764',
#         u'http://www.w3.org/2000/01/rdf-schema#label':
#           u'GIST with decreased sensitivity to therapy',
#         u'id':
#           u'http://ohsu.edu/cgd/87752f6c'},
#         u'phenotype_label':
#           u'GIST with decreased sensitivity to therapy',
#     u'sources':
#       u'http://www.ncbi.nlm.nih.gov/pubmed/18955458'
#  }]
        # annotation keys
        TYPE = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type'
        LABEL = 'http://www.w3.org/2000/01/rdf-schema#label'
        # useful
        # ECO_0000033 traceable author statement
        # RO_0002558 has evidence
        # RO_0002200 has phenotype
        # RO_0002606 is substance that treats
        # SO_0001059 sequence_alteration
        # BFO_0000159 has quality
        # OMIM_606764

        feature = association['feature']

        f = protocol.Feature()
        f.featureType = protocol.OntologyTerm.fromJsonDict({
            "term": feature[TYPE],
            "id": feature['id'],
            "sourceVersion": self._version,
            "sourceName": self._getPrefix(
                self._getPrefixURL(association['id']))
        })
        # TODO connect with real feature Ids
        f.id = feature['id']
        f.referenceName = feature[LABEL]
        f.attributes = protocol.Attributes.fromJsonDict(
            {"vals":  feature})
        f.childIds = []

        fpa = protocol.FeaturePhenotypeAssociation()
        fpa.id = association['id']
        fpa.features = [f]
        msg = 'Association: genotype:[{}] phenotype:[{}] environment:[{}] ' \
              'evidence:[{}] publications:[{}]'
        fpa.description = msg.format(
            association['feature_label'],
            association['phenotype_label'],
            association['environment_label'],
            self._getIdentifier(association['evidence']),
            association['sources']
            )
        evidence = protocol.Evidence()
        phenotype = association['phenotype']
        evidence.evidenceType = protocol.OntologyTerm.fromJsonDict({
            "term": association['evidence_type'],
            "id": phenotype['id'],
            "sourceVersion": self._version,
            "sourceName": self._getPrefix(
                self._getPrefixURL(association['id']))
        })
        evidence.description = self._getIdentifier(association['evidence'])
        # TODO there is nowhere in evidence to place list of sources?
        fpa.evidence = [evidence]

        # map environment (drug)
        environmentalContext = protocol.EnvironmentalContext()
        environment = association['environment']
        environmentalContext.id = environment['id']
        environmentalContext.description = association['environment_label']
        envType = protocol.OntologyTerm.fromJsonDict({
            "id": 'http://purl.obolibrary.org/obo/RO_0002606',
            "term": environment['id'],
            "sourceVersion": self._version,
            "sourceName": self._getPrefix(
                self._getPrefixURL(association['id']))
        })
        environmentalContext.environmentType = envType
        fpa.environmentalContexts = [environmentalContext]

        phenotypeInstance = protocol.PhenotypeInstance()
        phenotypeInstance.type = protocol.OntologyTerm.fromJsonDict({
            "term": phenotype[TYPE],
            "id": phenotype['id'],
            "sourceVersion": self._version,
            "sourceName": self._getPrefix(
                self._getPrefixURL(association['id']))
        })
        phenotypeInstance.description = phenotype[LABEL]
        fpa.phenotype = phenotypeInstance
        fpa.phenotypeAssociationSetId = self.getId()

# # DBG
# if not protocol.FeaturePhenotypeAssociation.validate(fpa.toJsonDict()):
#     # import pdb; pdb.set_trace() # DBG
#     import pprint
#     pp = pprint.PrettyPrinter(indent=2)
#     pp.pprint(f.toJsonString())
#     raise exceptions.RequestValidationFailureException(
#         fpa.toJsonDict(), protocol.FeaturePhenotypeAssociation)

        return fpa
