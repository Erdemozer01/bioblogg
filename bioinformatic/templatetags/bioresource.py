from django import template
from bioinformatic.models import BiologicalResourcesDatabases

register = template.Library()


@register.simple_tag()
def resource():
    return BiologicalResourcesDatabases.objects.all()

