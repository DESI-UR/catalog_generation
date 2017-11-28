
import catalog_generator

CatalogGenerator = catalog_generator.CatalogGenerator(configFile="data/catgen.cfg")
CatalogGenerator.getConfig()
CatalogGenerator.getLists()
CatalogGenerator.generate_randoms()
CatalogGenerator.generate_mocks()

