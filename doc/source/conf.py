import sphinx_rtd_theme


html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

def setup(app):
    app.add_css_file("main_stylesheet.css")

extensions = ['breathe']
breathe_projects = { 'analyticalscattering': '../xml' }
templates_path = ['_templates']
html_static_path = ['_static']
source_suffix = '.rst'
project = 'Analytical Scattering'
copyright = '2020, Ecole Polytechnique, Marc Bakry'
author = 'Marc Bakry'
# html_logo = 'castor_logo.png'

exclude_patterns = []