"""
Module responsible for translating g2p data into GA4GH native
objects.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import json
import rdflib

import ga4gh.datamodel as datamodel
import ga4gh.exceptions as exceptions
import ga4gh.protocol as protocol
import ga4gh.datamodel.sequenceAnnotations as sequenceAnnotations
import ga4gh.datamodel.genotype_phenotype as g2p
import ga4gh.pb as pb

# annotation keys
TYPE = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type'
LABEL = 'http://www.w3.org/2000/01/rdf-schema#label'
HAS_QUALITY = 'http://purl.obolibrary.org/obo/BFO_0000159'



class PhenotypeAssociationFeatureSet(g2p.G2PUtility,
                                     sequenceAnnotations.Gff3DbFeatureSet):
    """
    An rdf object store.  The cancer genome database
    [Clinical Genomics Knowledge Base]
    (http://nif-crawler.neuinfo.org/monarch/ttl/cgd.ttl),
    published by the Monarch project, was the source of Evidence.
    """

    def __init__(self, parentContainer, localId):
        super(PhenotypeAssociationFeatureSet, self).__init__(parentContainer, localId)

    # mimic featureset
    def populateFromRow(self, row):
        """
        Populates the instance variables of this FeatureSet from the specified
        DB row.
        """
        self._dbFilePath = row[b'dataUrl']
        self.populateFromFile(self._dbFilePath)

    def populateFromFile(self, dataUrl):
        """
        Populates the instance variables of this FeatureSet from the specified
        data URL.
        Initialize dataset, using the passed dict of sources
        [{source,format}] see rdflib.parse() for more
        If path is set, this backend will load itself
        """
        self._dbFilePath = dataUrl

        # initialize graph
        self._rdfGraph = rdflib.ConjunctiveGraph()
        # save the path
        self._dataUrl = dataUrl

        try:
            self._scanDataFiles(self._dataUrl, ['*.ttl', '*.xml'])
        except AttributeError:
            pass

        # extract version
        cgdTTL = rdflib.URIRef("http://data.monarchinitiative.org/ttl/cgd.ttl")
        versionInfo = rdflib.URIRef(
            u'http://www.w3.org/2002/07/owl#versionInfo')
        self._version = None
        for s, p, o in self._rdfGraph.triples((cgdTTL, versionInfo, None)):
            self._version = o.toPython()

    # mimic featureset
    def getFeature(self, compoundId):
        """
        find a feature and return ga4gh representation
        """
        print("getFeature.... {}".format(compoundId.featureId))
        featureRef = rdflib.URIRef(compoundId.featureId)
        featureDetails =  self._detailTuples([featureRef])
        feature = {}
        for f in featureDetails:
            feature[f['predicate']] = []

        for f in featureDetails:
            feature[f['predicate']].append(f['object'])

        f = protocol.Feature()

        term = protocol.OntologyTerm()
        # Schema for feature only supports one type of `type`
        # here we default to first OBO defined
        for featureType in feature[TYPE]:
            if "obolibrary" in featureType:
                term.term = featureType
                term.id = self._getPrefixURL(featureType)
                f.feature_type.MergeFrom(term)
                break ;


        f.id = str(compoundId)
        # Schema for feature only supports one type of `name` `symbol`
        # here we default to shortest for symbol and longest for name
        feature[LABEL].sort(key = len)
        f.gene_symbol = feature[LABEL][0]
        f.name = feature[LABEL][-1]

        f.attributes.MergeFrom(protocol.Attributes())
        for key in feature:
            for val in feature[key]:
                f.attributes.vals[key].values.add().string_value = val

        return f


    # mimic featureset
    def getFeatures(self, referenceName=None, start=None, end=None,
                    pageToken=None, pageSize=None,
                    featureTypes=None, parentId=None,
                    name=None, geneSymbol=None, numFeatures=10):

        # query to do search
        query = self._filterSearchFeaturesRequest(referenceName,geneSymbol,name)
        associations = self._rdfGraph.query(query)
        # associations is now a dict with rdflib terms with variable and
        # URIrefs or literals

        # given get the details for the feature,phenotype and environment
        associations_details = self._detailTuples(
                                    self._extractAssociationsDetails(
                                        associations))

        # association_details is now a list of {subject,predicate,object}
        # for each of the association detail
        # http://nmrml.org/cv/v1.0.rc1/doc/doc/objectproperties/BFO0000159___-324347567.html
        # label "has quality at all times" (en)
        associationList = []
        for assoc in associations.bindings:
            if '?feature' in assoc:
                association = self._bindingsToDict(assoc)
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

        # create GA4GH objects
        associations = [self._toGA4GH(assoc) for
                        assoc in associationList]
        features = []
        nextPageToken = None
        for annotation in associations:
            for feature in annotation.features:
                yield feature, (
                    str(nextPageToken)
                    if nextPageToken is not None else None)

    def _filterSearchFeaturesRequest(self, reference_name,gene_symbol,name):

        filters = []
        query = self._baseQuery()

        filters = []

        if reference_name:
            filters.append('regex(?feature_label, "{}")'
                           .format("NOT SUPPORTED"))

        if gene_symbol:
            filters.append('regex(?feature_label, "{}")'
                           .format(gene_symbol))
        if name:
            filters.append('regex(?feature_label, "{}")'
                           .format(name))

        # apply filters
        filter = "FILTER ({})".format(' && '.join(filters))
        if len(filters) == 0:
            filter = ""
        query = query.replace("#%FILTER%", filter)

        print(query) # TODO cleanup
        return query
